"""Measurements, local annotation glyphs, meshes, and particle models."""

import itertools
from copy import deepcopy

import numpy as np

from .data import finite
from .mesh import merge, unit
from .molecule import axes_of, ellipsoid
from .primitives import cylinder, dashed, indexed, sphere, sweep, text_mesh
from .scene import Geometry, PickTarget
from .themes import categorical, rgb, tables


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
    color = rgb(p.get("color", 0xFFA500 if rep == "shape-orientation" else 0x90EE90))
    radius = float(
        p.get("linesSize", p.get("sizeFactor", 0.075 if rep == "distance" else 0.04))
    )
    result = Geometry()
    atoms = targets(state, rep)
    pieces = []
    if rep in ("shape-orientation", "shape-plane"):
        center, axes, extent = axes_of(xyz)
        extent *= float(p.get("scaleFactor", 1))
        if rep == "shape-plane":
            vertices = (
                center
                + np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) * extent[:2] @ axes[:2]
            )
            pieces.append(indexed(vertices, [[0, 1, 2], [0, 2, 3]], color))
        else:
            visuals = p.get("visuals", ["box"])
            radius *= float(p.get("radiusScale", 2))
            if "axes" in visuals:
                pieces.extend(
                    cylinder(center - a * r, center + a * r, radius, color, 0)
                    for a, r in zip(axes, extent)
                )
            if "box" in visuals:
                vertices = (
                    center
                    + np.array(list(itertools.product((-1, 1), repeat=3)))
                    * extent
                    @ axes
                )
                pieces.extend(
                    cylinder(vertices[i], vertices[j], radius, color, 0)
                    for i in range(8)
                    for j in (i ^ 1, i ^ 2, i ^ 4)
                    if i < j
                )
            if "ellipsoid" in visuals:
                pieces.append(ellipsoid(center, axes, extent, color, 0, 16))
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
        defaults = {
            "distance": ["lines", "text"],
            "angle": ["vectors", "sector", "text"],
            "dihedral": ["extenders", "arms", "sector", "text"],
        }
        visuals = p.get("visuals", defaults[rep])
        if any(v in visuals for v in ("lines", "vectors", "arms")):
            pieces.extend(
                dashed(
                    a,
                    b,
                    radius,
                    color,
                    0,
                    count=max(
                        1,
                        int(
                            np.linalg.norm(b - a)
                            / (
                                2
                                * float(
                                    p.get(
                                        "dashLength", 0.2 if rep == "distance" else 0.04
                                    )
                                )
                            )
                        ),
                    ),
                )
                for a, b in itertools.pairwise(xyz)
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
                arm_length = min(
                    np.linalg.norm(xyz[0] - center), np.linalg.norm(xyz[2] - center)
                )
                if np.linalg.norm(np.cross(a, b)) < 1e-7:
                    axis = unit(np.cross(a, np.eye(3)[np.argmin(abs(a))]))
            else:
                axis = unit(xyz[2] - xyz[1])
                center = (xyz[1] + xyz[2]) / 2
                a = unit((xyz[0] - xyz[1]) - axis * np.dot(xyz[0] - xyz[1], axis))
                b = unit((xyz[3] - xyz[2]) - axis * np.dot(xyz[3] - xyz[2], axis))
                angle = np.arctan2(np.dot(np.cross(a, b), axis), np.dot(a, b))
                arm_length = min(
                    np.linalg.norm(np.cross(xyz[0] - xyz[1], axis)),
                    np.linalg.norm(np.cross(xyz[3] - xyz[2], axis)),
                )
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
            arc = center + arm_length * float(p.get("arcScale", 0.7)) * (
                np.cos(t)[:, None] * a + np.sin(t)[:, None] * np.cross(axis, a)
            )
            if "arc" in visuals:
                pieces.extend(
                    cylinder(a, b, radius, color, 0, 6)
                    for a, b in itertools.pairwise(arc)
                )
            if "sector" in visuals and abs(angle) > 1e-7:
                sector = indexed(
                    np.r_[center[None, :], arc],
                    [[0, i, i + 1] for i in range(1, len(arc))],
                    color,
                    0,
                    float(p.get("sectorOpacity", 0.75)),
                )
                result.add(sector, atoms, True)
            if rep == "dihedral" and "extenders" in visuals:
                pieces.extend(
                    dashed(a, b, radius, color, 0)
                    for a, b in ((xyz[0], arc[0]), (xyz[3], arc[-1]))
                )
            center = arc[len(arc) // 2]
        if p.get("label", True) and "text" in visuals:
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
    if (
        not np.isfinite(values).all()
        or min(values[:3]) <= 0
        or any(not 0 < a < 180 for a in values[3:])
    ):
        raise ValueError(
            "Unit cell must have positive lengths and nondegenerate angles"
        )
    cell = gemmi.UnitCell(*values)
    if not np.isfinite(cell.volume) or cell.volume <= 0:
        raise ValueError("Unit cell must have positive volume")
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
            vertex_rows = [r for r in tables["mesh_vertex"] if str(r["mesh_id"]) == key]
            lookup = {int(r["vertex_id"]): i for i, r in enumerate(vertex_rows)}
            if len(lookup) != len(vertex_rows):
                raise ValueError("Duplicate mesh vertex_id")
            vertices = [[float(r[k]) for k in ("x", "y", "z")] for r in vertex_rows]
            try:
                ids = [
                    lookup[int(r["vertex_id"])]
                    for r in tables["mesh_triangle"]
                    if str(r["mesh_id"]) == key
                ]
            except KeyError as exc:
                raise ValueError("Mesh triangle refers to a missing vertex_id") from exc
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


def particle_target(target, p, base=None):
    """Prepare reusable local target geometry before applying particle transforms."""
    from .data import Grid, local_path, read_grid
    from .volume import isosurface, molecular_grid, volume_geometry

    kind = target.get("kind", "shape")
    if kind == "shape":
        geometry = mesh_geometry(target, p)
    elif kind == "volume":
        grid = target.get("grid")
        if isinstance(grid, str):
            grid = read_grid(local_path(grid, base))
        elif isinstance(grid, dict):
            grid = Grid(grid["values"], grid.get("transform", np.eye(4)))
        if not isinstance(grid, Grid):
            raise ValueError("Particle volume target requires a local grid")
        geometry = volume_geometry(
            grid,
            target.get("type", "isosurface"),
            {**p, **target.get("params", {})},
            "medium",
        )
    elif kind == "structure":
        xyz = finite(target.get("positions"), (None, 3), "target atom positions")
        radii = np.broadcast_to(target.get("radii", 1), len(xyz))
        if not np.isfinite(radii).all() or (radii <= 0).any():
            raise ValueError("Target atom radii must be positive")
        color = rgb(target.get("color", p.get("color", 0xCCCCCC)))
        if target.get("type", "spacefill") == "blob-surface":
            grid, level = molecular_grid(xyz, radii, "gaussian-surface", p, "medium")
            mesh = isosurface(grid, level, color)
        else:
            mesh = merge(
                [sphere(pos, rad, color, 0, 12) for pos, rad in zip(xyz, radii)]
            )
        geometry = Geometry().add(mesh)
    else:
        raise ValueError("Particle target kind must be shape, structure, or volume")
    if not geometry.pieces or geometry.volumes:
        raise ValueError(
            "Particle target requires a nonempty mesh or dot representation"
        )
    vertices = np.concatenate([piece.mesh.vertices for piece in geometry.pieces])
    center = finite(
        target.get("center", (vertices.min(axis=0) + vertices.max(axis=0)) / 2),
        (3,),
        "target center",
    )
    return geometry, center, kind


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
    prepared = {}
    for i, (row, color) in enumerate(zip(rows, palette)):
        color = rgb(p.get("color", row.get("color", color)))
        radius = float(row.get("radius", 1)) * float(p.get("sizeFactor", 1))
        if not np.isfinite(radius) or radius <= 0:
            raise ValueError("Particle radius and sizeFactor must be positive")
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
            if not np.allclose(rotation @ rotation.T, np.eye(3), atol=1e-5):
                raise ValueError("Particle axes must be orthonormal")
        scale = finite(row.get("scale", [1, 1, 1]), (3,), "particle scale")
        if (scale <= 0).any():
            raise ValueError("Particle scale must be positive")
        pieces = []
        if rep == "particle-spacefill":
            pieces.append(
                ellipsoid(
                    center,
                    rotation,
                    scale * radius,
                    color,
                    0,
                    12,
                )
            )
        elif rep == "particle-orientation":
            pieces.extend(
                cylinder(
                    center,
                    center + a * float(p.get("axisLength", 10)),
                    0.02,
                    rgb(p.get(k, default)),
                    0,
                )
                for a, k, default in zip(
                    rotation,
                    ("xColor", "yColor", "zColor"),
                    (0xFF0000, 0x008000, 0x0000FF),
                )
            )
        elif rep == "particle-fibers":
            from scipy.interpolate import CubicHermiteSpline

            points = finite(row.get("points"), (None, 3), "fiber points")
            if len(points) < 2:
                raise ValueError("Fiber requires at least two points")
            segments = int(p.get("linearSegments", 8))
            if not 1 <= segments <= 128:
                raise ValueError("linearSegments must be 1..128")
            tangents = np.gradient(points, axis=0)
            tangents[[0, -1]] *= 0.5
            xyz = CubicHermiteSpline(np.arange(len(points)), points, tangents)(
                np.linspace(0, len(points) - 1, (len(points) - 1) * segments + 1)
            )
            widths = np.full(len(xyz), radius * float(p.get("tubeSizeFactor", 0.5)))
            pieces.append(
                sweep(
                    xyz,
                    widths,
                    widths,
                    np.zeros_like(xyz),
                    np.tile(color, (len(xyz), 1)),
                    np.zeros(len(xyz), int),
                    "ellipse",
                    12,
                )
            )
        else:
            target_id = row.get("target")
            if not isinstance(target_id, (str, int)) or str(target_id) not in data.get(
                "targets", {}
            ):
                raise ValueError(
                    "Particle target must name an entry in targets; XYZ arrow endpoints are not targets"
                )
            key = str(target_id)
            if key not in prepared:
                prepared[key] = particle_target(
                    data["targets"][key], p, data.get("_base")
                )
            geometry, origin, kind = prepared[key]
            multiplier = radius if p.get("scaleByRadius", kind == "shape") else 1
            transform = np.diag(scale * multiplier) @ rotation
            for piece in geometry.pieces:
                mesh = deepcopy(piece.mesh)
                mesh.vertices = ((mesh.vertices - origin) @ transform + center).astype(
                    np.float32
                )
                mesh.normals = unit(mesh.normals @ np.linalg.inv(transform).T).astype(
                    np.float32
                )
                if np.linalg.det(transform) < 0:
                    mesh.faces = mesh.faces[:, ::-1].copy()
                mesh.owners[:] = 0
                if p.get("targetColor", "particle") == "particle":
                    mesh.colors[:] = color
                pieces.append(mesh)
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
        for row in steps:
            colors = [
                rgb(
                    p.get(
                        "color",
                        row.get(
                            key + "Color",
                            tables()["color_maps"]["NtCColors"].get(
                                str(row["class"]) + "_" + suffix, 0xFFA10A
                            ),
                        ),
                    )
                )
                for key, suffix in (("upper", "Upr"), ("lower", "Lwr"))
            ]
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
                verts[1] = (verts[0] + verts[4]) / 2
                for faces, color in zip(
                    (
                        [[0, 2, 3], [0, 3, 1], [0, 1, 2]],
                        [[4, 3, 2], [4, 1, 3], [4, 2, 1]],
                    ),
                    colors,
                ):
                    pieces.append(
                        indexed(
                            verts[np.array(faces)].reshape(-1, 3),
                            np.arange(9).reshape(-1, 3),
                            color,
                            0,
                        )
                    )
            else:
                pieces.extend(
                    cylinder(
                        a,
                        b,
                        float(p.get("sizeFactor", 0.35)),
                        colors[int(k >= (len(xyz) - 1) / 2)],
                        0,
                    )
                    for k, (a, b) in enumerate(itertools.pairwise(xyz))
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
