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
    if np.linalg.det(transform) < 0:
        m.faces = m.faces[:, ::-1].copy()
    return m


def axes_of(points):
    center = points.mean(axis=0)
    _, _, axes = np.linalg.svd(points - center, full_matrices=True)
    if np.linalg.det(axes) < 0:
        axes[-1] *= -1
    projected = (points - center) @ axes.T
    lo, hi = projected.min(axis=0), projected.max(axis=0)
    center = center + (lo + hi) / 2 @ axes
    extent = np.maximum((hi - lo) / 2, 0.2)
    return center, axes, extent


def atomic(state, representation, c, r, p, quality):
    detail = QUALITIES[quality][1]
    factor = float(
        p.get(
            "sizeFactor",
            {"spacefill": 1, "ellipsoid": 1, "point": 1, "line": 2}.get(
                representation, 0.15
            ),
        )
    )
    if not np.isfinite(factor) or factor <= 0:
        raise ValueError("sizeFactor must be positive")
    visible = np.array(
        [
            not p.get("ignoreHydrogens", False) or a.element.upper() not in ("H", "D")
            for a in state.atoms
        ]
    )
    if p.get("traceOnly", False):
        visible &= np.array(
            [
                a.name == "CA"
                if a.kind == "protein"
                else a.name in ("C4'", "P")
                if a.kind == "nucleic"
                else False
                for a in state.atoms
            ]
        )
    stride = int(p.get("stride", 1))
    if stride < 1:
        raise ValueError("stride must be positive")
    if representation == "point":
        visible &= np.arange(len(state.atoms)) % stride == 0
    defaults = {
        "line": ["intra-bond", "inter-bond", "element-point", "element-cross"],
        "point": ["element-point"],
        "ellipsoid": ["ellipsoid-mesh", "intra-bond", "inter-bond"],
        "ball-and-stick": ["element-sphere", "intra-bond", "inter-bond"],
        "spacefill": ["element-sphere"],
    }
    visuals = [
        v.removeprefix("structure-")
        for v in p.get("visuals", defaults.get(representation, []))
    ]
    if representation == "ellipsoid" and not any(
        len(a.aniso) == 6 and any(a.aniso) for a in state.atoms
    ):
        raise ValueError("Ellipsoid requires anisotropic displacement tensors")
    bonded = {i for bond in state.bonds for i in bond}
    neighbors = {i: [] for i in bonded}
    for i, j in state.bonds:
        neighbors[i].append(j)
        neighbors[j].append(i)
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
            if "element-cross" in visuals and (
                p.get("crosses", "lone") == "all" or i not in bonded
            ):
                half_cross = float(p.get("crossSize", 0.35)) / 2
                for axis in np.eye(3):
                    pieces.append(
                        cylinder(
                            pos - axis * half_cross,
                            pos + axis * half_cross,
                            r[i] * factor * 0.02,
                            c[i],
                            i,
                            6,
                        )
                    )
            if representation == "ellipsoid":
                if len(a.aniso) != 6 or not any(a.aniso):
                    continue
                u11, u22, u33, u12, u13, u23 = a.aniso
                values, axes = np.linalg.eigh(
                    [[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]]
                )
                # Match the pinned Mol* multiplier; an explicit probability is an extension.
                multiplier = 1.5958
                if "probability" in p:
                    from scipy.stats import chi2

                    probability = float(p["probability"])
                    if not 0 < probability < 1:
                        raise ValueError("probability must be between zero and one")
                    multiplier = np.sqrt(chi2.ppf(probability, 3))
                pieces.append(
                    ellipsoid(
                        pos,
                        axes.T,
                        np.sqrt(np.abs(values)) * multiplier * factor,
                        c[i],
                        i,
                        detail,
                    )
                )
            elif representation == "polyhedron":
                pass
            elif "element-point" in visuals or "element-sphere" in visuals:
                pieces.append(
                    sphere(
                        pos,
                        radius * 0.075 if "element-point" in visuals else radius,
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
                min(r[i], r[j]) * factor * 0.02
                if representation == "line"
                else min(r[i], r[j])
                * factor
                * float(
                    p.get(
                        "sizeAspectRatio",
                        0.1 if representation == "ellipsoid" else 2 / 3,
                    )
                )
            )
            order = int(getattr(state.model.bond[k], "order", 1))
            multiple = p.get(
                "multipleBonds", "offset" if representation == "line" else "symmetric"
            )
            count = order if order in (2, 3) and multiple != "off" else 1
            normal = unit(
                np.cross(
                    unit(b - a), [1, 0, 0] if abs(unit(b - a)[0]) < 0.8 else [0, 1, 0]
                )
            )
            for neighbor in neighbors.get(i, []) + neighbors.get(j, []):
                if neighbor in (i, j) or state.atoms[neighbor].element.upper() in (
                    "H",
                    "D",
                ):
                    continue
                reference = state.coords[neighbor] - a
                projected = reference - unit(b - a) * np.dot(reference, unit(b - a))
                if np.linalg.norm(projected) > 1e-6:
                    normal = unit(projected)
                    break
            link_scale = float(
                p.get("linkScale", 0.5 if representation == "line" else 0.45)
            )
            spacing = float(
                p.get("linkSpacing", 0.1 if representation == "line" else 1)
            )
            multi_radius = radius * link_scale / (0.5 * count)
            for v in range(count):
                if count == 1:
                    shift, bond_radius = 0, radius
                elif multiple == "offset":
                    shift = (0, 1, -1)[v] * (
                        radius + multi_radius + link_scale * radius * spacing
                    )
                    bond_radius = radius if v == 0 else multi_radius
                else:
                    shift = (
                        (v - (count - 1) / 2)
                        * (2 if count == 2 else 1)
                        * (radius - multi_radius)
                        * spacing
                    )
                    bond_radius = multi_radius
                offset = normal * shift
                mid = (a + b) / 2 + offset
                for start, end, color, owner in (
                    (a + offset, mid, c[i], i),
                    (mid, b + offset, c[j], j),
                ):
                    pieces.append(
                        cylinder(
                            start,
                            end,
                            bond_radius,
                            color,
                            owner,
                            detail,
                        )
                    )
            if order == 4 and p.get("aromaticBonds", True):
                aromatic_scale = float(p.get("aromaticScale", 0.3))
                offset = (
                    normal
                    * radius
                    * (1 + aromatic_scale * (1 + float(p.get("aromaticSpacing", 1.5))))
                )
                pieces.append(
                    dashed(
                        a + offset,
                        b + offset,
                        radius * aromatic_scale,
                        c[i],
                        i,
                        count=int(p.get("aromaticDashCount", 2)),
                    )
                )
    return merge(pieces)


def polymer(state, representation, c, r, p, quality):
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
            representation,
            [
                "polymer-trace",
                "polymer-gap",
                "nucleotide-ring",
                "nucleotide-atomic-ring-fill",
                "nucleotide-atomic-bond",
                "nucleotide-atomic-element",
            ],
        ),
    )
    group_lookup = {i: g for g in groups for i in g}
    size = float(p.get("sizeFactor", 0.3 if representation == "backbone" else 0.2))
    if not np.isfinite(size) or size <= 0:
        raise ValueError("sizeFactor must be positive")
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
            result.append(
                sphere(
                    state.coords[ids[0]], size * r[ids[0]], c[ids[0]], ids[0], radial
                )
            )
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
                        cylinder(state.coords[i], middle, size * r[i], c[i], i, radial),
                        cylinder(middle, state.coords[j], size * r[j], c[j], j, radial),
                    ]
                )
            if "polymer-backbone-sphere" in visuals:
                result.extend(
                    sphere(state.coords[i], size * r[i], c[i], i, radial) for i in ids
                )
            continue
        t = np.arange(len(ids))
        tangent = np.gradient(points, axis=0)
        samples = np.linspace(0, len(ids) - 1, (len(ids) - 1) * linear + 1)
        xyz = CubicHermiteSpline(t, points, tangent)(samples)
        owner = np.array(ids)[
            np.clip(np.floor(samples + 0.5).astype(int), 0, len(ids) - 1)
        ]
        widths = size * np.interp(samples, t, r[ids])
        thickness = widths.copy()
        hints = np.zeros_like(xyz)
        for k, i in enumerate(owner):
            a = state.atoms[i]
            if representation == "putty" and "bfactorScale" in p:
                widths[k] = thickness[k] = size * max(
                    0.25, np.sqrt(max(a.bfactor, 0) / float(p.get("bfactorScale", 25)))
                )
            elif representation != "putty" and (
                a.ss in ("H", "S") or a.kind == "nucleic"
            ):
                widths[k] *= float(p.get("aspectRatio", 5))
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
    for group, names in (
        nucleotide if any(v.startswith("nucleotide-") for v in visuals) else []
    ):
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
    theme = p.get(
        "sizeTheme",
        "uncertainty"
        if representation == "putty"
        else "uniform"
        if representation
        in ("ellipsoid", "cartoon", "backbone", "point", "line", "orientation")
        else "physical",
    )
    r = sizes(state.atoms, theme, p.get("sizeParams", {}), data)
    result = Geometry()
    if representation in ("cartoon", "backbone", "putty"):
        return result.add(polymer(state, representation, c, r, p, quality), state.atoms)
    if representation == "carbohydrate":
        from .carbohydrate import geometry

        return geometry(state, p, quality, c if mode != "auto" else None)
    if representation in ("molecular-surface", "gaussian-surface", "gaussian-volume"):
        from .scene import Volume
        from .volume import isosurface, molecular_grid, transfer_function, wireframe

        grid, level = molecular_grid(state.coords, r, representation, p, quality)
        if representation == "gaussian-volume":
            color_grid = np.empty((*grid.values.shape, 3), np.float32)
            tree = cKDTree(state.coords)
            yz = np.indices(grid.values.shape[1:]).reshape(2, -1).T
            for x in range(grid.values.shape[0]):
                points = grid.world(np.c_[np.full(len(yz), x), yz])
                color_grid[x] = c[tree.query(points)[1]].reshape(
                    *grid.values.shape[1:], 3
                )
            result.volumes.append(
                Volume(
                    grid,
                    transfer_function(grid, p),
                    float(p.get("alpha", 1)),
                    min(grid.spacing),
                    color_grid=color_grid,
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
    if representation == "plane":
        from .structure_plane import plane

        return result.add(plane(state, c, r, p), state.atoms, True)
    if representation == "orientation":
        center, axes, extent = axes_of(state.coords)
        extent *= float(p.get("sizeFactor", 1))
        visuals = p.get("visuals", ["orientation-ellipsoid-mesh"])
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
