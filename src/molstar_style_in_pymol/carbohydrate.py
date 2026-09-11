"""SNFG symbols and covalent glycan links, using Mol* residue tables."""

import numpy as np

from .mesh import merge
from .primitives import convex, cylinder, indexed, sphere
from .scene import Geometry
from .themes import rgb, tables


def sugar_info(resn):
    item = tables()["saccharides"].get(resn, {})
    return item.get("shape", "FlatHexagon"), item.get("color", 0xF1ECE1)


def symbol(center, shape, radius, color, owner, axes=None, detail=12):
    axes = np.eye(3) if axes is None else np.asarray(axes)
    if shape == "FilledSphere":
        return sphere(center, radius, color, owner, detail)
    if "Cube" in shape or shape == "FlatBox":
        v = (
            np.array(
                [[x, y, z] for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)], float
            )
            * 0.65
        )
        if shape == "FlatBox":
            v[:, 2] *= 0.35
    elif "Diamond" in shape and "Prism" not in shape:
        v = np.r_[np.eye(3), -np.eye(3)]
        if shape == "FlatDiamond":
            v[:, 2] *= 0.35
    elif "Cone" in shape:
        return cylinder(
            center - radius * axes[1],
            center + radius * axes[1],
            radius,
            color,
            owner,
            detail,
            0,
        )
    else:
        count = {
            "Pentagon": 5,
            "PentagonalPrism": 5,
            "DiamondPrism": 4,
            "HexagonalPrism": 6,
            "HeptagonalPrism": 7,
            "FilledStar": 10,
        }.get(shape, 6)
        t = np.arange(count) * 2 * np.pi / count
        r = (
            np.where(np.arange(count) % 2, 0.45, 1)
            if shape == "FilledStar"
            else np.ones(count)
        )
        ring = np.c_[r * np.cos(t), r * np.sin(t)]
        v = np.array([[x, y, z] for z in (-0.2, 0.2) for x, y in ring])
        if shape == "FilledStar":
            verts = np.r_[v, [[0, 0, -0.2], [0, 0, 0.2]]]
            faces = []
            for i in range(count):
                j = (i + 1) % count
                faces += [
                    [2 * count, j, i],
                    [2 * count + 1, i + count, j + count],
                    [i, j, j + count],
                    [i, j + count, i + count],
                ]
            return indexed(center + radius * verts @ axes, faces, color, owner)
    mesh = convex(center + radius * v @ axes, color, owner)
    if shape in ("CrossedCube", "DividedDiamond"):
        relative = (mesh.vertices - center) @ axes.T
        mask = (
            (relative[:, 0] * relative[:, 1] > 0)
            if shape == "CrossedCube"
            else relative[:, 1] > 0
        )
        mesh.colors[mask] = rgb(0xF1ECE1)
    return mesh


def geometry(state, params, quality, colors=None):
    groups = {}
    for i, a in enumerate(state.atoms):
        if a.resn in tables()["saccharides"]:
            groups.setdefault((a.model, a.segi, a.chain, a.resi, a.resn), []).append(i)
    result = Geometry()
    if not groups:
        raise ValueError("No recognized carbohydrate residues in the selection")
    centers = {}
    pieces = []
    visuals = params.get(
        "visuals",
        ["carbohydrate-symbol", "carbohydrate-link", "carbohydrate-terminal-link"],
    )
    radius = float(params.get("sizeFactor", 1.75))
    for key, ids in groups.items():
        xyz = state.coords[ids]
        center = xyz.mean(axis=0)
        _, _, axes = np.linalg.svd(xyz - center, full_matrices=True)
        shape, color = sugar_info(key[-1])
        if np.linalg.det(axes) < 0:
            axes[-1] *= -1
        if "carbohydrate-symbol" in visuals:
            pieces.append(
                symbol(
                    center,
                    shape,
                    radius,
                    rgb(color) if colors is None else colors[ids[0]],
                    ids[0],
                    axes,
                )
            )
        for i in ids:
            centers[i] = (center, ids[0])
    links = set()
    for i, j in state.bonds:
        if i not in centers and j not in centers:
            continue
        visual = (
            "carbohydrate-link"
            if i in centers and j in centers
            else "carbohydrate-terminal-link"
        )
        if visual not in visuals:
            continue
        a, ai = centers.get(i, (state.coords[i], i))
        b, bi = centers.get(j, (state.coords[j], j))
        key = tuple(sorted((ai, bi)))
        if ai != bi and key not in links:
            links.add(key)
            pieces.append(
                cylinder(
                    a, b, float(params.get("linkSizeFactor", 0.25)), rgb(0x999999), ai
                )
            )
    return result.add(merge(pieces), state.atoms)
