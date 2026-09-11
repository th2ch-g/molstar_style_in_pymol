"""Scalar-field surfaces, slices, segmentations, and true density sampling."""

import numpy as np
from scipy.ndimage import distance_transform_edt, gaussian_filter
from scipy.spatial import cKDTree

from .data import Grid
from .mesh import Mesh, merge, unit
from .primitives import cylinder, sphere
from .scene import Geometry, PickTarget, Volume
from .themes import categorical, rgb, scale


def grid_box(coords, radii, spacing, budget=512 * 1024**2):
    if not np.isfinite(spacing) or spacing <= 0:
        raise ValueError("resolution must be finite and positive")
    low = np.min(coords - radii[:, None], axis=0) - spacing * 2
    high = np.max(coords + radii[:, None], axis=0) + spacing * 2
    shape = np.maximum(3, np.ceil((high - low) / spacing).astype(int) + 1)
    if int(np.prod(shape)) * 24 > budget:
        raise ValueError(
            "Surface grid exceeds cache_mb; increase resolution spacing or reduce selection"
        )
    transform = np.eye(4)
    transform[:3, :3] *= spacing
    transform[:3, 3] = low
    return Grid(np.zeros(shape, np.float32), transform)


def molecular_grid(coords, radii, kind, params, quality):
    spacing = float(
        params.get(
            "resolution", {"low": 1.0, "medium": 0.7, "high": 0.5}.get(quality, 0.7)
        )
    )
    radii = np.asarray(radii, float) + float(params.get("radiusOffset", 0))
    smooth = float(params.get("smoothness", 1.5))
    probe = float(params.get("probeRadius", 1.4))
    if np.any(radii <= 0) or smooth <= 0 or probe < 0:
        raise ValueError(
            "Radii and smoothness must be positive; probeRadius must be nonnegative"
        )
    inflated = radii + probe if kind == "molecular-surface" else radii * 2
    grid = grid_box(
        coords, inflated, spacing, int(params.get("_budget", 512 * 1024**2))
    )
    origin = grid.transform[:3, 3]
    for center, radius, bound in zip(coords, radii, inflated):
        lo = np.maximum(0, np.floor((center - bound - origin) / spacing).astype(int))
        hi = np.minimum(
            grid.values.shape,
            np.ceil((center + bound - origin) / spacing).astype(int) + 1,
        )
        points = np.meshgrid(
            *(np.arange(a, b) * spacing + o for a, b, o in zip(lo, hi, origin)),
            indexing="ij",
            sparse=True,
        )
        d2 = sum((p - c) ** 2 for p, c in zip(points, center))
        region = tuple(slice(a, b) for a, b in zip(lo, hi))
        if kind == "molecular-surface":
            grid.values[region] = np.maximum(grid.values[region], d2 <= bound**2)
        else:
            grid.values[region] += np.exp(-smooth * d2 / radius**2).astype(np.float32)
    if kind == "molecular-surface":
        # Euclidean erosion of the solvent-accessible solid retains reentrant pockets.
        grid.values = distance_transform_edt(grid.values > 0, sampling=spacing).astype(
            np.float32
        )
        level = max(probe, spacing * 0.5)
    else:
        level = np.exp(-smooth)
    return grid, level


def isosurface(grid, level, color=(0.2, 0.6, 0.8), coords=None, colors=None):
    from skimage.measure import marching_cubes

    level = grid.level(level)
    if not float(grid.values.min()) < level < float(grid.values.max()):
        raise ValueError(
            f"Isovalue {level:g} must lie inside [{grid.values.min():g}, {grid.values.max():g}]"
        )
    vertices, faces, normals, _ = marching_cubes(
        grid.values, level, allow_degenerate=False
    )
    vertices = grid.world(vertices)
    normals = unit(normals @ np.linalg.inv(grid.transform[:3, :3]))
    if np.linalg.det(grid.transform[:3, :3]) < 0:
        faces = faces[:, ::-1]
    if coords is not None and len(coords):
        _, owners = cKDTree(coords).query(vertices)
        vertex_colors = np.asarray(colors)[owners]
    else:
        owners = np.zeros(len(vertices), int)
        vertex_colors = np.tile(rgb(color), (len(vertices), 1))
    return Mesh(vertices, normals, vertex_colors, faces, owners)


def wireframe(mesh, radius=0.035, detail=6):
    edges = np.unique(
        np.sort(
            np.concatenate(
                [mesh.faces[:, [0, 1]], mesh.faces[:, [1, 2]], mesh.faces[:, [2, 0]]]
            ),
            axis=1,
        ),
        axis=0,
    )
    return merge(
        [
            cylinder(
                mesh.vertices[a],
                mesh.vertices[b],
                radius,
                mesh.colors[a],
                int(mesh.owners[a]),
                detail,
            )
            for a, b in edges
        ]
    )


def transfer_function(grid, params):
    value = params.get("transferFunction")
    if value is not None:
        rows = np.array([[float(r[0]), *rgb(r[1]), float(r[2])] for r in value])
        if (
            len(rows) < 2
            or not np.isfinite(rows).all()
            or np.any(np.diff(rows[:, 0]) <= 0)
            or np.any((rows[:, 4] < 0) | (rows[:, 4] > 1))
        ):
            raise ValueError(
                "transferFunction requires increasing [value, color, opacity] rows"
            )
        return rows
    lo, hi = float(grid.values.min()), float(grid.values.max())
    control = np.asarray(
        params.get(
            "controlPoints", [[0, 0], [0.25, 0], [0.5, 0.08], [0.75, 0.3], [1, 0.8]]
        ),
        float,
    )
    if (
        control.ndim != 2
        or control.shape[1] != 2
        or len(control) < 2
        or np.any(np.diff(control[:, 0]) <= 0)
        or np.any((control < 0) | (control > 1))
    ):
        raise ValueError(
            "controlPoints requires increasing normalized [value, opacity] pairs"
        )
    colors = scale(control[:, 0], params.get("colorList", "viridis"), [0, 1])
    return np.c_[lo + control[:, 0] * (hi - lo), colors, control[:, 1]]


def transfer_values(volume, values, step=1):
    t = volume.transfer
    out = np.stack(
        [
            np.interp(
                values,
                t[:, 0],
                t[:, i],
                left=0 if i == 4 else t[0, i],
                right=0 if i == 4 else t[-1, i],
            )
            for i in range(1, 5)
        ],
        axis=-1,
    )
    out[..., 3] = (
        1 - np.power(1 - out[..., 3], step / max(volume.step, 1e-6))
    ) * volume.opacity
    return out


def slice_mesh(grid, params):
    axis = "xyz".index(params.get("dimension", "z"))
    other = [i for i in range(3) if i != axis]
    index = float(
        params.get(
            "index", params.get("absoluteIndex", (grid.values.shape[axis] - 1) / 2)
        )
    )
    if "relativeIndex" in params:
        index = float(params["relativeIndex"]) * (grid.values.shape[axis] - 1)
    if not 0 <= index <= grid.values.shape[axis] - 1:
        raise ValueError("Slice index is outside the grid")
    axes = [np.arange(grid.values.shape[i]) for i in other]
    a, b = np.meshgrid(*axes, indexing="ij")
    positions = np.zeros((a.size, 3))
    positions[:, axis] = index
    positions[:, other[0]], positions[:, other[1]] = a.ravel(), b.ravel()
    vertices = grid.world(positions)
    c = scale(
        grid.sample(vertices), params.get("colorList", "viridis"), params.get("domain")
    )
    ni, nj = a.shape
    ids = np.arange(ni * nj).reshape(ni, nj)
    x, y, z, w = (
        ids[:-1, :-1].ravel(),
        ids[1:, :-1].ravel(),
        ids[:-1, 1:].ravel(),
        ids[1:, 1:].ravel(),
    )
    faces = np.concatenate([np.c_[x, y, z], np.c_[y, w, z]])
    normal = unit(np.cross(grid.transform[:3, other[0]], grid.transform[:3, other[1]]))
    return Mesh(
        vertices,
        np.tile(normal, (len(vertices), 1)),
        c,
        faces,
        np.zeros(len(vertices), int),
    )


def volume_geometry(grid, representation, params, quality="medium"):
    geometry = Geometry()
    target = (PickTarget(label=grid.label),)
    color = rgb(params.get("color", 0x33AADD))
    if representation == "isosurface":
        iso = params.get("isoValue", {"kind": "relative", "relativeValue": 1})
        mesh = isosurface(grid, iso, color)
        visuals = params.get("visuals", ["solid"])
        if params.get("showFaces", True) and (
            "solid" in visuals or "isosurface" in visuals
        ):
            geometry.add(mesh, target)
        if params.get("showWireframe") or "wireframe" in visuals:
            geometry.add(
                wireframe(mesh, float(params.get("sizeFactor", 0.035))), target
            )
    elif representation == "slice":
        geometry.add(slice_mesh(grid, params), target, True)
    elif representation == "dot":
        level = grid.level(
            params.get("isoValue", {"kind": "relative", "relativeValue": 1})
        )
        indices = np.argwhere(grid.values >= level)
        stride = max(1, int(params.get("stride", 1)))
        indices = indices[::stride]
        if len(indices) > int(params.get("maxPoints", 100_000)):
            raise ValueError("Too many volume dots; increase stride or isovalue")
        vertices = grid.world(indices)
        colors = scale(
            grid.values[tuple(indices.T)], params.get("colorList", "viridis")
        )
        radius = float(params.get("sizeFactor", 0.12))
        geometry.add(
            merge([sphere(p, radius, c, 0, 6) for p, c in zip(vertices, colors)]),
            target,
        )
    elif representation == "segment":
        if not np.allclose(grid.values, np.rint(grid.values)):
            raise ValueError("Segmentation requires an integer label grid")
        ids = params.get("segments", [int(v) for v in np.unique(grid.values) if v != 0])
        colors = categorical(ids, params.get("colorList", "many-distinct"))
        for value, color in zip(ids, colors):
            mask = grid.values == value
            if not mask.any():
                raise ValueError(f"Segment {value} is absent from the label grid")
            # Padding closes segments touching an outer face of the input grid.
            transform = grid.transform.copy()
            transform[:3, 3] -= transform[:3, :3] @ np.ones(3)
            field = gaussian_filter(
                np.pad(mask.astype(np.float32), 1), float(params.get("smoothness", 0))
            )
            mesh = isosurface(Grid(field, transform), 0.5, color)
            geometry.add(mesh, (PickTarget(label=f"{grid.label}: segment {value}"),))
    elif representation in ("direct-volume", "gaussian-volume"):
        geometry.volumes.append(
            Volume(
                grid,
                transfer_function(grid, params),
                float(params.get("alpha", 1)),
                float(params.get("step", grid.spacing.min() * 0.5)),
            )
        )
    else:
        raise ValueError(f"Unknown volume representation {representation}")
    return geometry


def density_slices(volume, quality="medium", direction=None):
    """Sample transfer-function density into transparent quads for native ray.

    A dedicated export supplies a camera direction. Retained proxies use three
    orthogonal stacks with one-third optical depth, avoiding a disappearing slab
    when the camera is parallel to one stack.
    """
    from .primitives import indexed

    grid = volume.grid
    if direction is None:
        copies = []
        for axis in np.eye(3):
            piece = density_slices(volume, quality, axis)
            piece.alphas = 1 - (1 - piece.alphas) ** (1 / 3)
            copies.append(piece)
        return merge(copies)
    normal = unit(direction)
    x = unit(np.cross(normal, np.eye(3)[np.argmin(abs(normal))]))
    y = np.cross(normal, x)
    basis = np.stack([x, y, normal])
    corners = (
        np.array(np.meshgrid(*[[0, n - 1] for n in grid.values.shape], indexing="ij"))
        .reshape(3, -1)
        .T
    )
    projected = grid.world(corners) @ basis.T
    low, high = projected.min(axis=0), projected.max(axis=0)
    samples = {"low": 12, "medium": 20, "high": 32, "highest": 48}.get(quality, 20)
    spacing = max(
        grid.spacing.min() * (0.25 if quality in ("high", "highest") else 1),
        float(np.max(high - low)) / samples,
    )
    dims = np.maximum(2, np.ceil((high - low) / spacing).astype(int) + 1)
    axes = [np.linspace(a, b, n) for a, b, n in zip(low, high, dims)]
    xx, yy = np.meshgrid(axes[0], axes[1], indexing="ij")
    vertices, colors, alphas, faces = [], [], [], []
    for z in axes[2]:
        positions = np.c_[xx.ravel(), yy.ravel(), np.full(xx.size, z)] @ basis
        rgba = transfer_values(volume, grid.sample(positions), spacing)
        ids = np.arange(xx.size).reshape(xx.shape)
        a, b, c, d = (
            ids[:-1, :-1].ravel(),
            ids[1:, :-1].ravel(),
            ids[:-1, 1:].ravel(),
            ids[1:, 1:].ravel(),
        )
        triangles = np.concatenate([np.c_[a, b, c], np.c_[b, d, c]])
        keep = rgba[triangles, 3].max(axis=1) > 0.001
        triangles = triangles[keep]
        if not len(triangles):
            continue
        offset = len(vertices)
        vertices.extend(positions)
        colors.extend(rgba[:, :3])
        alphas.extend(rgba[:, 3])
        faces.extend(triangles + offset)
    if not faces:
        return Mesh([], [], [], [], [])
    mesh = indexed(vertices, faces, colors, 0)
    mesh.alphas = np.asarray(alphas, np.float32)
    return mesh


def density_projection(volume):
    """Three pre-integrated density planes for PyMOL's single-surface ray mode."""
    meshes = []
    for axis in range(3):
        field = np.moveaxis(volume.grid.values, axis, 0)
        total = np.zeros((*field.shape[1:], 4), float)
        for values in field:
            rgba = transfer_values(volume, values, float(volume.grid.spacing[axis]))
            weight = (1 - total[:, :, 3]) * rgba[:, :, 3]
            total[:, :, :3] += rgba[:, :, :3] * weight[:, :, None]
            total[:, :, 3] += weight
        alpha = total[:, :, 3].ravel()
        colors = np.divide(
            total[:, :, :3].reshape(-1, 3),
            alpha[:, None],
            out=np.zeros((alpha.size, 3)),
            where=alpha[:, None] > 1e-8,
        )
        mesh = slice_mesh(volume.grid, {"dimension": "xyz"[axis]})
        mesh.colors = colors.astype(np.float32)
        mesh.alphas = (1 - (1 - alpha) ** (1 / 3)).astype(np.float32)
        mesh.faces = mesh.faces[mesh.alphas[mesh.faces].max(axis=1) > 0.005]
        meshes.append(mesh)
    return merge(meshes)
