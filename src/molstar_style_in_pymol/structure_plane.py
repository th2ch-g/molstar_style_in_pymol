"""Sample an atom-colored cross section, independently of fitted shape planes."""

import numpy as np
from scipy.spatial.transform import Rotation

from .data import finite
from .mesh import Mesh, unit
from .themes import rgb


def plane(state, colors, radii, params):
    from .molecule import axes_of

    center, axes, extent = axes_of(state.coords)
    if params.get("frame", "principalAxes") == "boundingBox":
        lo, hi = state.coords.min(axis=0), state.coords.max(axis=0)
        center, axes, extent = (lo + hi) / 2, np.eye(3), (hi - lo) / 2
    axis = "abc".index(params.get("axis", "c"))
    other = [i for i in range(3) if i != axis]
    major, minor, normal = axes[other[0]], axes[other[1]], axes[axis]
    half = extent[other].copy()
    depth = extent[axis]
    rotation = params.get("rotation", {})
    angle = float(rotation.get("angle", 0))
    sphere_radius = np.linalg.norm(state.coords - center, axis=1).max()
    if params.get("extent") == "sphere" or angle:
        half[:] = np.sqrt(3) * sphere_radius
        depth = sphere_radius
    margin = float(params.get("margin", 4))
    half += margin
    depth += margin
    if angle:
        ra = unit(finite(rotation.get("axis", [1, 0, 0]), (3,), "rotation axis") @ axes)
        if np.linalg.norm(ra) < 0.9:
            raise ValueError("Plane rotation axis cannot be zero")
        major, minor, normal = Rotation.from_rotvec(ra * np.deg2rad(angle)).apply(
            [major, minor, normal]
        )
    if params.get("mode", "frame") == "plane":
        definition = params.get("plane", {})
        center = finite(definition.get("point", [0, 0, 0]), (3,), "plane point")
        normal = unit(finite(definition.get("normal", [1, 0, 0]), (3,), "plane normal"))
        if np.linalg.norm(normal) < 0.9:
            raise ValueError("Plane normal cannot be zero")
        major = unit(np.cross(normal, np.eye(3)[np.argmin(abs(normal))]))
        minor = np.cross(normal, major)
    center = center + normal * depth * float(params.get("offset", 0))
    resolution = float(params.get("imageResolution", 0.5))
    if not np.isfinite(resolution) or resolution <= 0 or np.any(half <= 0):
        raise ValueError("Plane resolution and extent must be positive")
    shape = np.maximum(1, np.ceil(2 * half / resolution).astype(int))
    if np.prod(shape) * 256 > int(params.get("_budget", 512 * 1024**2)):
        raise ValueError("Plane image exceeds cache_mb; increase imageResolution")
    spacing = 2 * half / shape
    x, y = np.meshgrid(
        *(
            np.linspace(-h + d / 2, h - d / 2, n)
            for h, d, n in zip(half, spacing, shape)
        ),
        indexing="ij",
    )
    positions = center + x.ravel()[:, None] * major + y.ravel()[:, None] * minor
    background = rgb(params.get("defaultColor", 0xCCCCCC))
    pixel_colors = np.tile(background, (len(positions), 1))
    pixel_alpha = np.full(len(positions), 0.0 if params.get("cutout", False) else 1.0)
    owners = np.zeros(len(positions), int)
    nearest = np.full(len(positions), np.inf)
    for i, (pos, radius) in enumerate(zip(state.coords, radii)):
        distance = abs(np.dot(pos - center, normal))
        if distance > radius:
            continue
        delta = positions - pos
        distances = np.linalg.norm(delta, axis=1)
        tolerance = (
            resolution * np.cos(distance / radius)
            if params.get("antialias", True)
            else 0
        )
        keep = (distances <= radius + tolerance / 2) & (distances < nearest)
        coverage = (
            np.clip((radius - distances[keep]) / max(tolerance, 1e-12) + 0.5, 0, 1)
            if tolerance
            else np.ones(keep.sum())
        )
        pixel_colors[keep] = (
            colors[i]
            if params.get("cutout", False)
            else background + coverage[:, None] * (colors[i] - background)
        )
        pixel_alpha[keep] = coverage if params.get("cutout", False) else 1
        nearest[keep], owners[keep] = distances[keep], i
    keep = pixel_alpha > 0
    if params.get("extent") == "sphere":
        keep &= np.linalg.norm(positions - center, axis=1) <= sphere_radius + margin
    positions = positions[keep]
    corners = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) * spacing / 2
    vertices = (
        positions[:, None, :] + corners[:, :1] * major + corners[:, 1:] * minor
    ).reshape(-1, 3)
    faces = (
        np.arange(len(positions))[:, None, None] * 4 + np.array([[0, 1, 2], [0, 2, 3]])
    ).reshape(-1, 3)
    return Mesh(
        vertices,
        np.tile(unit(np.cross(major, minor)), (len(vertices), 1)),
        np.repeat(pixel_colors[keep], 4, axis=0),
        faces,
        np.repeat(owners[keep], 4),
        alphas=np.repeat(pixel_alpha[keep], 4),
    )
