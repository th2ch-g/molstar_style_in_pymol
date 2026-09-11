"""Measurements, local annotation glyphs, meshes, and particle models."""

import itertools

import numpy as np

from .data import finite
from .mesh import merge, unit
from .molecule import axes_of, ellipsoid
from .primitives import arrow, cylinder, dashed, indexed, sphere, text_mesh
from .scene import Geometry, PickTarget
from .themes import categorical, rgb


def targets(state, label):
    return state.atoms if state and state.atoms else (PickTarget(label=label),)


def positions(state, data, count=None):
    values = data.get("positions")
    if values is None and state:
        indices = data.get("indices", list(range(len(state.atoms))))
        if any(int(i) != i or not 0 <= i < len(state.atoms) for i in indices):
            raise ValueError("indices must refer to zero-based selected atoms")
        values = state.coords[indices]
    values = finite(values, (count, 3), "positions")
    return values


def measurements(state, rep, p, data):
    n = {"distance": 2, "angle": 3, "dihedral": 4}.get(rep)
    xyz = positions(state, data, n)
    color = rgb(p.get("color", 0xFFFF00))
    radius = float(p.get("sizeFactor", 0.06))
    result = Geometry()
    atoms = targets(state, rep)
    pieces = []
    if rep in ("shape-orientation", "shape-plane"):
        center, axes, extent = axes_of(xyz)
        if rep == "shape-plane":
            vertices = (
                center
                + np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) * extent[:2] @ axes[:2]
            )
            pieces.append(indexed(vertices, [[0, 1, 2], [0, 2, 3]], color))
        else:
            pieces.extend(
                arrow(center, center + a * r, radius, c, 0)
                for a, r, c in zip(axes, extent, np.eye(3))
            )
    elif rep in ("shape-label", "custom-label", "annotation-label"):
        text = data.get("text", p.get("text"))
        if text is None:
            raise ValueError("Label requires text")
        pieces.append(
            text_mesh(
                str(text),
                xyz.mean(axis=0),
                color,
                float(p.get("textSize", 0.6)),
                0,
                font=p.get("font"),
            )
        )
    else:
        pieces.extend(
            dashed(a, b, radius, color, 0) for a, b in itertools.pairwise(xyz)
        )
        if rep == "distance":
            value = np.linalg.norm(xyz[1] - xyz[0])
            label = f"{value:.2f} A"
            center = xyz.mean(axis=0)
        else:
            if rep == "angle":
                a, b = unit(xyz[0] - xyz[1]), unit(xyz[2] - xyz[1])
                center = xyz[1]
                angle = np.arccos(np.clip(np.dot(a, b), -1, 1))
                axis = unit(np.cross(a, b))
            else:
                axis = unit(xyz[2] - xyz[1])
                center = (xyz[1] + xyz[2]) / 2
                a = unit((xyz[0] - xyz[1]) - axis * np.dot(xyz[0] - xyz[1], axis))
                b = unit((xyz[3] - xyz[2]) - axis * np.dot(xyz[3] - xyz[2], axis))
                angle = np.arctan2(np.dot(np.cross(a, b), axis), np.dot(a, b))
            value = np.degrees(angle)
            label = f"{value:.1f} deg"
            if (
                np.linalg.norm(axis) < 1e-7
                or np.linalg.norm(a) < 1e-7
                or np.linalg.norm(b) < 1e-7
            ):
                raise ValueError(
                    "Measurement is undefined for coincident or collinear points"
                )
            t = np.linspace(0, angle, 33)
            arc = center + float(p.get("arcScale", 1)) * (
                np.cos(t)[:, None] * a + np.sin(t)[:, None] * np.cross(axis, a)
            )
            pieces.extend(
                cylinder(a, b, radius, color, 0, 6) for a, b in itertools.pairwise(arc)
            )
            center = arc[len(arc) // 2]
        if p.get("label", True):
            pieces.append(
                text_mesh(
                    data.get("text", label),
                    center,
                    color,
                    float(p.get("textSize", 0.5)),
                    0,
                )
            )
        result.metadata["value"] = float(value)
    return result.add(merge(pieces), atoms, True)


def unitcell(state, p, data):
    import gemmi

    values = data.get("cell")
    if values is None or len(values) != 6:
        raise ValueError("unitcell requires cell: [a,b,c,alpha,beta,gamma]")
    cell = gemmi.UnitCell(*values)
    basis = np.array(cell.orth.mat)
    origin = np.asarray(data.get("origin", [0, 0, 0]))
    vertices = (
        origin
        + np.array([[x, y, z] for x in (0, 1) for y in (0, 1) for z in (0, 1)])
        @ basis.T
    )
    pieces = [
        cylinder(
            vertices[i],
            vertices[j],
            float(p.get("sizeFactor", 0.08)),
            rgb(p.get("color", 0xDD8800)),
            0,
        )
        for i in range(8)
        for j in (i ^ 1, i ^ 2, i ^ 4)
        if i < j
    ]
    return Geometry().add(merge(pieces), targets(state, "unit cell"))


def mesh_geometry(data, p):
    if "tables" in data:
        tables = data["tables"]
        rows = []
        for group in tables.get("mesh", []):
            key = str(group["id"])
            vertices = [
                [float(r[k]) for k in ("x", "y", "z")]
                for r in tables["mesh_vertex"]
                if str(r["mesh_id"]) == key
            ]
            ids = [
                int(r["vertex_id"])
                for r in tables["mesh_triangle"]
                if str(r["mesh_id"]) == key
            ]
            rows.append(
                {
                    "vertices": vertices,
                    "faces": np.array(ids).reshape(-1, 3),
                    "label": key,
                }
            )
        data = {**data, "meshes": rows}
    rows = data.get("meshes", [data])
    result = Geometry()
    for i, row in enumerate(rows):
        vertices = finite(row.get("vertices"), (None, 3), "vertices")
        faces = np.asarray(row.get("faces"))
        if (
            faces.ndim != 2
            or faces.shape[1] != 3
            or not np.isfinite(faces).all()
            or np.any(faces != np.rint(faces))
        ):
            raise ValueError("faces must contain integer triangles")
        transform = finite(row.get("transform", np.eye(4)), (4, 4), "transform")
        vertices = vertices @ transform[:3, :3].T + transform[:3, 3]
        colors = (
            [rgb(v) for v in row["colors"]]
            if "colors" in row
            else rgb(row.get("color", p.get("color", 0x33AADD)))
        )
        mesh = indexed(
            vertices, faces.astype(int), colors, 0, float(row.get("opacity", 1))
        )
        result.add(mesh, (PickTarget(label=str(row.get("label", i))),))
    return result


def particles(data, rep, p):
    result = Geometry()
    rows = data.get("particles")
    if rows is None:
        coords = finite(data.get("positions"), (None, 3), "positions")
        radii = np.broadcast_to(data.get("radii", 1), len(coords))
        rows = [{"position": v, "radius": r} for v, r in zip(coords, radii)]
    palette = categorical(
        [r.get(p.get("colorBy", "entity"), i) for i, r in enumerate(rows)],
        p.get("colorList", "many-distinct"),
    )
    if p.get("scalarColor"):
        from .themes import scale

        palette = scale(
            [r[p.get("colorBy", "value")] for r in rows],
            p.get("colorList", "viridis"),
            p.get("colorParams", {}).get("domain"),
        )
    for i, (row, color) in enumerate(zip(rows, palette)):
        color = rgb(p.get("color", row.get("color", color)))
        radius = float(row.get("radius", 1)) * float(p.get("sizeFactor", 1))
        center = finite(row.get("position", [0, 0, 0]), (3,), "particle position")
        rotation = np.eye(3)
        if "quaternion" in row:
            from scipy.spatial.transform import Rotation

            rotation = (
                Rotation.from_quat(
                    finite(row["quaternion"], (4,), "quaternion [x,y,z,w]")
                )
                .as_matrix()
                .T
            )
        if "axes" in row:
            rotation = finite(row["axes"], (3, 3), "axes")
        pieces = []
        if rep == "particle-spacefill":
            pieces.append(
                ellipsoid(
                    center,
                    rotation,
                    np.array(row.get("scale", [radius] * 3)),
                    color,
                    0,
                    12,
                )
            )
        elif rep == "particle-orientation":
            pieces.extend(
                arrow(center, center + a * radius, 0.08 * radius, c, 0)
                for a, c in zip(rotation, np.eye(3))
            )
        elif rep == "particle-fibers":
            points = finite(row.get("points"), (None, 3), "fiber points")
            pieces.extend(
                cylinder(a, b, radius, color, 0) for a, b in itertools.pairwise(points)
            )
            pieces.extend(sphere(a, radius, color, 0, 8) for a in points)
        else:
            target = finite(row.get("target"), (3,), "particle target")
            pieces.append(arrow(center, target, radius * 0.1, color, 0))
            pieces.append(sphere(target, radius * 0.25, color, 0, 8))
        result.add(merge(pieces), (PickTarget(label=str(row.get("label", i))),))
    return result


def annotated(state, rep, p, data):
    result = Geometry()
    atoms = targets(state, rep)
    pieces = []
    if rep == "membrane-orientation":
        membrane = data.get("membrane")
        if membrane is None:
            raise ValueError(
                "membrane-orientation requires membrane.center, normal, thickness, radius"
            )
        center = finite(membrane["center"], (3,), "membrane center")
        normal = unit(finite(membrane["normal"], (3,), "membrane normal"))
        if np.linalg.norm(normal) < 0.9:
            raise ValueError("Membrane normal cannot be zero")
        a = unit(np.cross(normal, [1, 0, 0] if abs(normal[0]) < 0.8 else [0, 1, 0]))
        b = np.cross(normal, a)
        t = np.linspace(0, 2 * np.pi, 65)
        ring = (np.cos(t)[:, None] * a + np.sin(t)[:, None] * b) * float(
            membrane["radius"]
        )
        for sign in (-1, 1):
            pos = center + normal * float(membrane["thickness"]) / 2 * sign
            verts = np.r_[pos[None, :], pos + ring]
            mesh = indexed(
                verts,
                [[0, i, i + 1] for i in range(1, len(ring))],
                rgb(p.get("color", 0xAAAAAA)),
                0,
                float(p.get("alpha", 0.35)),
            )
            result.add(mesh, atoms)
            pieces.extend(
                cylinder(a + pos, b + pos, 0.08, rgb(0xCCCC33), 0)
                for a, b in itertools.pairwise(ring)
            )
    elif rep == "assembly-symmetry":
        symmetry = data.get("symmetry")
        if symmetry is None:
            raise ValueError(
                "assembly-symmetry requires local symmetry.axes and optional cage"
            )
        for ax in symmetry.get("axes", []):
            color = rgb(ax.get("color", 0xFF8833))
            start = finite(ax["start"], (3,), "axis start")
            end = finite(ax["end"], (3,), "axis end")
            pieces.append(
                cylinder(start, end, float(p.get("sizeFactor", 0.2)), color, 0)
            )
            pieces.append(text_mesh(str(ax.get("order", 2)), end, color, 0.6, 0))
        if "cage" in symmetry:
            vertices = finite(symmetry["cage"]["vertices"], (None, 3), "cage vertices")
            for i, j in symmetry["cage"]["edges"]:
                pieces.append(
                    cylinder(vertices[i], vertices[j], 0.08, rgb(0xAAAAAA), 0)
                )
    elif rep in ("confal-pyramids", "ntc-tube"):
        steps = data.get("steps")
        if not steps:
            raise ValueError(
                "DNATCO display requires local steps with positions, class, and score"
            )
        palette = categorical([r["class"] for r in steps])
        for row, color in zip(steps, palette):
            xyz = finite(row["positions"], (None, 3), "step positions")
            if rep == "confal-pyramids":
                if len(xyz) != 5:
                    raise ValueError(
                        "Confal positions must be O3-prime, P, OP1, OP2, O5-prime"
                    )
                scale = (float(row["score"]) - 20) / 100
                verts = xyz.copy()
                verts[[0, 2, 3, 4]] = xyz[1] + (xyz[[0, 2, 3, 4]] - xyz[1]) * (
                    1 + scale
                )
                pieces.append(
                    indexed(
                        verts,
                        [
                            [0, 2, 3],
                            [0, 3, 1],
                            [0, 1, 2],
                            [4, 3, 2],
                            [4, 1, 3],
                            [4, 2, 1],
                        ],
                        color,
                        0,
                    )
                )
            else:
                pieces.extend(
                    cylinder(a, b, float(p.get("sizeFactor", 0.35)), color, 0)
                    for a, b in itertools.pairwise(xyz)
                )
    elif rep == "tunnel":
        tunnels = data.get("tunnels", [data])
        for row in tunnels:
            points = finite(row.get("positions"), (None, 3), "tunnel positions")
            radii = finite(row.get("radii"), (len(points),), "tunnel radii")
            if (radii <= 0).any():
                raise ValueError("Tunnel radii must be positive")
            color = rgb(row.get("color", 0xFFAA33))
            pieces.extend(sphere(a, r, color, 0, 12) for a, r in zip(points, radii))
            pieces.extend(
                cylinder(a, b, r, color, 0, 12, s)
                for a, b, r, s in zip(points[:-1], points[1:], radii[:-1], radii[1:])
            )
    else:
        raise ValueError(f"Unknown annotation representation {rep}")
    return result.add(merge(pieces), atoms)
