"""Independent atomic, polymer, and coarse molecular tessellation."""

import itertools

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.spatial import cKDTree

from .mesh import merge, unit
from .primitives import (
    arrow,
    convex,
    cylinder,
    dashed,
    indexed,
    sphere,
    sweep,
    text_mesh,
)
from .registry import QUALITIES
from .scene import Geometry
from .themes import colors, sizes


def residues(state):
    groups = {}
    for i, a in enumerate(state.atoms):
        groups.setdefault((a.model, a.segi, a.chain, a.resi), []).append(i)
    return list(groups.values())


def ellipsoid(center, axes, radii, color, owner, detail):
    m = sphere([0, 0, 0], 1, color, owner, detail)
    transform = np.diag(np.maximum(radii, 1e-4)) @ axes
    m.vertices = (m.vertices @ transform + center).astype(np.float32)
    m.normals = unit(m.normals @ np.linalg.inv(transform).T).astype(np.float32)
    return m


def axes_of(points):
    center = points.mean(axis=0)
    _, _, axes = np.linalg.svd(points - center, full_matrices=True)
    extent = np.maximum(np.max(np.abs((points - center) @ axes.T), axis=0), 0.2)
    if np.linalg.det(axes) < 0:
        axes[-1] *= -1
    return center, axes, extent


def atomic(state, representation, c, r, p, quality):
    detail = QUALITIES[quality][1]
    factor = float(
        p.get("sizeFactor", 1 if representation in ("spacefill", "ellipsoid") else 0.15)
    )
    if not np.isfinite(factor) or factor <= 0:
        raise ValueError("sizeFactor must be positive")
    visible = np.array(
        [
            not p.get("ignoreHydrogens", False) or a.element.upper() not in ("H", "D")
            for a in state.atoms
        ]
    )
    visuals = [v.removeprefix("structure-") for v in p.get("visuals", [])]
    pieces = []
    draw_atoms = (
        representation not in ("line", "polyhedron") or p.get("pointStyle") == "circle"
    )
    if visuals:
        draw_atoms = any(
            v in visuals
            for v in (
                "element-sphere",
                "element-point",
                "element-cross",
                "ellipsoid-mesh",
                "ellipsoid",
            )
        )
    if draw_atoms:
        for i, (a, pos) in enumerate(zip(state.atoms, state.coords)):
            if not visible[i]:
                continue
            radius = r[i] * factor
            if "element-cross" in visuals:
                for axis in np.eye(3):
                    pieces.append(
                        cylinder(
                            pos - axis * factor, pos + axis * factor, 0.025, c[i], i, 6
                        )
                    )
            elif representation == "ellipsoid" and len(a.aniso) == 6 and any(a.aniso):
                u11, u22, u33, u12, u13, u23 = a.aniso
                values, axes = np.linalg.eigh(
                    [[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]]
                )
                if np.any(values <= 0):
                    raise ValueError(
                        f"Nonpositive anisotropic tensor at atom {a.index}"
                    )
                from scipy.stats import chi2

                probability = float(p.get("probability", 0.5))
                if not 0 < probability < 1:
                    raise ValueError("probability must be between zero and one")
                pieces.append(
                    ellipsoid(
                        pos,
                        axes.T,
                        np.sqrt(values * chi2.ppf(probability, 3)) * factor,
                        c[i],
                        i,
                        detail,
                    )
                )
            elif representation == "polyhedron":
                pass
            else:
                pieces.append(
                    sphere(
                        pos,
                        radius if representation != "point" else factor,
                        c[i],
                        i,
                        detail,
                    )
                )
    if representation == "polyhedron":
        metals = {"ZN", "FE", "MG", "MN", "CA", "CU", "CO", "NI", "CD", "NA", "K"}
        tree = cKDTree(state.coords)
        for i, a in enumerate(state.atoms):
            if a.element.upper() not in metals:
                continue
            ids = [
                j
                for j in tree.query_ball_point(
                    state.coords[i], float(p.get("maxDistance", 3))
                )
                if j != i and state.atoms[j].element.upper() in ("O", "N", "S", "CL")
            ]
            if len(ids) >= 4:
                pieces.append(convex(state.coords[ids], c[i], i))
        if not pieces:
            raise ValueError(
                "Polyhedron requires a metal with at least four coordinating atoms"
            )
    if representation in ("ball-and-stick", "line", "ellipsoid") and (
        not visuals or any("bond" in v for v in visuals)
    ):
        for k, (i, j) in enumerate(state.bonds):
            if not visible[i] or not visible[j]:
                continue
            ai, aj = state.atoms[i], state.atoms[j]
            bond_kind = (
                "intra-bond"
                if (ai.model, ai.chain, ai.segi) == (aj.model, aj.chain, aj.segi)
                else "inter-bond"
            )
            if visuals and bond_kind not in visuals:
                continue
            a, b = state.coords[i], state.coords[j]
            radius = (
                float(p.get("lineSizeAttenuation", 0.04))
                if representation == "line"
                else min(r[i], r[j]) * factor * float(p.get("sizeAspectRatio", 2 / 3))
            )
            order = int(getattr(state.model.bond[k], "order", 1))
            multiple = p.get("multipleBonds", "symmetric")
            count = order if order in (2, 3) and multiple != "off" else 1
            normal = unit(
                np.cross(
                    unit(b - a), [1, 0, 0] if abs(unit(b - a)[0]) < 0.8 else [0, 1, 0]
                )
            )
            for v in range(count):
                offset = (
                    normal
                    * (v if multiple == "offset" else v - (count - 1) / 2)
                    * radius
                    * 2.6
                )
                mid = (a + b) / 2 + offset
                for start, end, color, owner in (
                    (a + offset, mid, c[i], i),
                    (mid, b + offset, c[j], j),
                ):
                    pieces.append(
                        cylinder(
                            start,
                            end,
                            radius if count == 1 else radius * 0.65,
                            color,
                            owner,
                            detail,
                        )
                    )
            if order == 4 and p.get("aromaticBonds", True):
                pieces.append(
                    dashed(
                        a + normal * radius * 2,
                        b + normal * radius * 2,
                        radius * 0.45,
                        c[i],
                        i,
                    )
                )
    return merge(pieces)


def polymer(state, representation, c, p, quality):
    linear, radial = QUALITIES[quality]
    linear = int(p.get("linearSegments", linear))
    radial = int(p.get("radialSegments", radial))
    if not 1 <= linear <= 128 or not 3 <= radial <= 128:
        raise ValueError("linearSegments must be 1..128; radialSegments must be 3..128")
    groups = residues(state)
    anchors = []
    nucleotide = []
    for group in groups:
        byname = {
            state.atoms[i].name: i for i in group if state.atoms[i].alt in ("", "A")
        }
        anchor = (
            byname.get("CA")
            if state.atoms[group[0]].kind == "protein"
            else byname.get("C4'", byname.get("P"))
            if state.atoms[group[0]].kind == "nucleic"
            else None
        )
        if anchor is not None:
            anchors.append(anchor)
        if state.atoms[group[0]].kind == "nucleic":
            nucleotide.append((group, byname))
    chunks = []
    gaps = []
    for i in anchors:
        a = state.atoms[i]
        if chunks:
            j = chunks[-1][-1]
            b = state.atoms[j]
            same = (a.model, a.chain, a.segi) == (b.model, b.chain, b.segi)
            distance = np.linalg.norm(state.coords[i] - state.coords[j])
            if same and distance < (10 if a.kind == "nucleic" else 4.5):
                chunks[-1].append(i)
                continue
            if same:
                gaps.append((j, i))
        chunks.append([i])
    result = []
    defaults = {
        "backbone": [
            "polymer-backbone-cylinder",
            "polymer-backbone-sphere",
            "polymer-gap",
        ],
        "putty": ["polymer-tube", "polymer-gap"],
    }
    visuals = p.get(
        "visuals",
        defaults.get(
            representation, ["polymer-trace", "polymer-gap", "nucleotide-ring"]
        ),
    )
    group_lookup = {i: g for g in groups for i in g}
    size = float(p.get("sizeFactor", 0.2))
    for ids in (
        chunks
        if any(
            v in visuals
            for v in (
                "polymer-trace",
                "polymer-tube",
                "polymer-backbone-cylinder",
                "polymer-backbone-sphere",
            )
        )
        else []
    ):
        if len(ids) == 1:
            result.append(sphere(state.coords[ids[0]], size, c[ids[0]], ids[0], radial))
            continue
        points = state.coords[ids]
        if representation == "backbone":
            for i, j in (
                itertools.pairwise(ids)
                if "polymer-backbone-cylinder" in visuals
                else []
            ):
                middle = (state.coords[i] + state.coords[j]) / 2
                result.extend(
                    [
                        cylinder(state.coords[i], middle, size, c[i], i, radial),
                        cylinder(middle, state.coords[j], size, c[j], j, radial),
                    ]
                )
            if "polymer-backbone-sphere" in visuals:
                result.extend(
                    sphere(state.coords[i], size, c[i], i, radial) for i in ids
                )
            continue
        t = np.arange(len(ids))
        tangent = np.gradient(points, axis=0)
        samples = np.linspace(0, len(ids) - 1, (len(ids) - 1) * linear + 1)
        xyz = CubicHermiteSpline(t, points, tangent)(samples)
        owner = np.array(ids)[
            np.clip(np.floor(samples + 0.5).astype(int), 0, len(ids) - 1)
        ]
        widths = np.full(len(xyz), size)
        thickness = widths.copy()
        hints = np.zeros_like(xyz)
        for k, i in enumerate(owner):
            a = state.atoms[i]
            if representation == "putty":
                widths[k] = thickness[k] = size * max(
                    0.25, np.sqrt(max(a.bfactor, 0) / float(p.get("bfactorScale", 25)))
                )
            elif a.ss in ("H", "S") or a.kind == "nucleic":
                widths[k] = size * float(p.get("aspectRatio", 5))
                if a.ss == "H" and p.get("tubularHelices", False):
                    thickness[k] = widths[k] * 1.5
                    widths[k] = thickness[k]
                if a.ss == "S":
                    segment = int(np.clip(np.floor(samples[k]), 0, len(ids) - 2))
                    if (
                        state.atoms[ids[segment + 1]].ss != "S"
                        or segment == len(ids) - 2
                    ):
                        f = samples[k] - segment
                        widths[k] *= float(p.get("arrowFactor", 1.5)) * max(0.02, 1 - f)
            group = group_lookup[i]
            names = {state.atoms[j].name: j for j in group}
            if "O" in names and "C" in names:
                hints[k] = state.coords[names["O"]] - state.coords[names["C"]]
        # Parallel transport handles missing peptide-plane directions.
        result.append(
            sweep(
                xyz,
                widths,
                thickness,
                hints,
                c[owner],
                owner,
                "rectangle"
                if all(state.atoms[i].ss == "S" for i in ids)
                else "ellipse",
                radial,
            )
        )
    if "polymer-gap" in visuals:
        result.extend(
            dashed(state.coords[i], state.coords[j], size * 0.7, c[i], i)
            for i, j in gaps
        )
    for group, names in nucleotide:
        ids = [
            names[n]
            for n in ("N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4")
            if n in names
        ]
        if len(ids) < 3:
            continue
        center, axes, extent = axes_of(state.coords[ids])
        if "nucleotide-block" in visuals:
            verts = np.array(
                [[x, y, z] for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
            ) * np.maximum(extent, [1, 0.6, 0.15])
            result.append(convex(center + verts @ axes, c[ids[0]], ids[0]))
        elif any(
            v in visuals for v in ("nucleotide-ring", "nucleotide-atomic-ring-fill")
        ):
            planar = (state.coords[ids] - center) @ axes[:2].T
            from scipy.spatial import ConvexHull

            hull = ConvexHull(planar).vertices
            ring = state.coords[np.array(ids)[hull]]
            verts = np.r_[ring + axes[2] * 0.12, ring - axes[2] * 0.12]
            result.append(convex(verts, c[ids[0]], ids[0]))
        if any(
            v in visuals
            for v in ("nucleotide-atomic-bond", "nucleotide-atomic-element")
        ):
            keep = set(ids)
            if "nucleotide-atomic-element" in visuals:
                result.extend(
                    sphere(state.coords[i], size, c[i], i, radial) for i in ids
                )
            if "nucleotide-atomic-bond" in visuals:
                result.extend(
                    cylinder(
                        state.coords[i], state.coords[j], size * 0.5, c[i], i, radial
                    )
                    for i, j in state.bonds
                    if i in keep and j in keep
                )
        anchor = names.get("C4'", names.get("P"))
        if anchor is not None:
            result.append(
                cylinder(state.coords[anchor], center, size, c[anchor], anchor, radial)
            )
    if "direction-wedge" in visuals:
        for ids in chunks:
            if len(ids) > 1:
                result.append(
                    arrow(
                        state.coords[ids[-2]],
                        state.coords[ids[-1]],
                        size * 2,
                        c[ids[-1]],
                        ids[-1],
                    )
                )
    return merge(result)


def blob(state, c, p, quality):
    count = max(1, int(p.get("blobSize", 30)))
    n = max(1, int(np.ceil(len(state.atoms) / count)))
    if p.get("method", "grid") == "clustering":
        from scipy.cluster.vq import kmeans2

        _, labels = kmeans2(
            state.coords, n, iter=int(p.get("clusterIterations", 2)), minit="++", seed=0
        )
    else:
        step = max(1, np.ptp(state.coords, axis=0).max() / max(1, n ** (1 / 3)))
        _, labels = np.unique(
            np.floor((state.coords - state.coords.min(axis=0)) / step),
            axis=0,
            return_inverse=True,
        )
    meshes = []
    for label in np.unique(labels):
        ids = np.flatnonzero(labels == label)
        center, axes, extent = axes_of(state.coords[ids])
        extent += np.mean([state.atoms[i].vdw for i in ids]) + float(
            p.get("radiusOffset", 0)
        )
        if p.get("shape", "ellipsoid") == "spherical-harmonics" and len(ids) >= 6:
            from scipy.special import sph_harm_y

            v = state.coords[ids] - center
            rad = np.maximum(np.linalg.norm(v, axis=1), 1e-6)

            def basis(vectors, degree):
                rr = np.maximum(np.linalg.norm(vectors, axis=1), 1e-6)
                theta = np.arccos(np.clip(vectors[:, 2] / rr, -1, 1))
                phi = np.arctan2(vectors[:, 1], vectors[:, 0])
                return np.array(
                    [
                        sph_harm_y(degree_index, m, theta, phi).real
                        if m >= 0
                        else sph_harm_y(degree_index, -m, theta, phi).imag
                        for degree_index in range(degree + 1)
                        for m in range(-degree_index, degree_index + 1)
                    ]
                ).T

            degree = int(p.get("degree", 2))
            b = basis(v, degree)
            coefficients = np.linalg.solve(
                b.T @ b + np.eye(b.shape[1]) * float(p.get("regularization", 0.05)),
                b.T @ (rad + 1.5),
            )
            mesh = sphere([0, 0, 0], 1, c[ids[0]], int(ids[0]), QUALITIES[quality][1])
            radius = np.clip(
                basis(mesh.vertices, degree) @ coefficients,
                0.5,
                (rad.max() + 1.5) * 1.25,
            )
            meshes.append(
                indexed(
                    center + mesh.vertices * radius[:, None],
                    mesh.faces,
                    c[ids[0]],
                    int(ids[0]),
                )
            )
        else:
            meshes.append(
                ellipsoid(
                    center, axes, extent, c[ids[0]], int(ids[0]), QUALITIES[quality][1]
                )
            )
    return merge(meshes)


def geometry(state, representation, mode, p, data, quality):
    c = colors(
        state.atoms, state.coords, mode, representation, p.get("colorParams", {}), data
    )
    theme = p.get("sizeTheme", "physical")
    r = sizes(state.atoms, theme, p.get("sizeParams", {}), data)
    result = Geometry()
    if representation in ("cartoon", "backbone", "putty"):
        return result.add(polymer(state, representation, c, p, quality), state.atoms)
    if representation == "carbohydrate":
        from .carbohydrate import geometry

        return geometry(state, p, quality)
    if representation in ("molecular-surface", "gaussian-surface", "gaussian-volume"):
        from .scene import Volume
        from .volume import isosurface, molecular_grid, transfer_function, wireframe

        grid, level = molecular_grid(state.coords, r, representation, p, quality)
        if representation == "gaussian-volume":
            result.volumes.append(
                Volume(
                    grid,
                    transfer_function(grid, p),
                    float(p.get("alpha", 1)),
                    min(grid.spacing),
                )
            )
        else:
            m = isosurface(
                grid, p.get("isoValue", level), coords=state.coords, colors=c
            )
            visuals = p.get("visuals", ["mesh"])
            if not all("wireframe" in v for v in visuals):
                result.add(m, state.atoms)
            if any("wireframe" in v for v in visuals):
                result.add(wireframe(m), state.atoms)
        return result
    if representation == "blob-surface":
        from .volume import wireframe

        m = blob(state, c, p, quality)
        visuals = p.get("visuals", ["blob-surface-mesh"])
        if any("mesh" in v for v in visuals):
            result.add(m, state.atoms)
        if any("wireframe" in v for v in visuals):
            result.add(wireframe(m), state.atoms)
        return result
    if representation in ("orientation", "plane"):
        center, axes, extent = axes_of(state.coords)
        if representation == "plane":
            verts = (
                center
                + np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]])
                * (extent[:2] + float(p.get("margin", 0)))
                @ axes[:2]
            )
            return result.add(
                indexed(verts, [[0, 1, 2], [0, 2, 3]], c[0], 0), state.atoms
            )
        visuals = p.get("visuals", ["axes", "box", "ellipsoid"])
        if "ellipsoid" in visuals or "orientation-ellipsoid-mesh" in visuals:
            result.add(
                ellipsoid(center, axes, extent, c[0], 0, QUALITIES[quality][1]),
                state.atoms,
            )
        if "axes" in visuals:
            for ax, length, col in zip(axes, extent, np.eye(3)):
                result.add(
                    arrow(center, center + ax * length, 0.12, col, 0), state.atoms
                )
        if "box" in visuals:
            verts = (
                center
                + np.array(
                    [[x, y, z] for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
                )
                * extent
                @ axes
            )
            for i in range(8):
                for j in (i ^ 1, i ^ 2, i ^ 4):
                    if j > i:
                        result.add(
                            cylinder(verts[i], verts[j], 0.05, c[0], 0), state.atoms
                        )
        return result
    if representation == "label":
        level = p.get("level", "residue")
        groups = (
            [[i] for i in range(len(state.atoms))]
            if level == "element"
            else residues(state)
        )
        if level == "chain":
            chains = {}
            for group in groups:
                chains.setdefault(state.atoms[group[0]].chain, []).extend(group)
            groups = list(chains.values())
        for group in groups:
            i = group[0]
            a = state.atoms[i]
            text = (
                a.name
                if level == "element"
                else a.chain
                if level == "chain"
                else f"{a.resn} {a.resi}"
            )
            result.add(
                text_mesh(
                    text,
                    state.coords[group].mean(axis=0),
                    c[i],
                    float(p.get("sizeFactor", 0.6)),
                    i,
                    font=p.get("font"),
                ),
                state.atoms,
                True,
            )
        return result
    return result.add(
        atomic(state, representation, c, r, p, quality),
        state.atoms,
        representation in ("point", "line"),
    )
