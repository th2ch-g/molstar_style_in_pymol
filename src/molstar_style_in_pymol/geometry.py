"""Representation dispatch and Mol*-inspired composition presets."""

import numpy as np

from . import analysis, formats, molecule, orbitals, panels, shapes, volume
from .registry import PARTICLE, PRESETS, SHAPE, STRUCTURE, VOLUME
from .scene import Geometry
from .themes import tables


def subset_data(data, mask):
    return {
        **data,
        "atom_data": {
            k: np.asarray(v)[mask] for k, v in data.get("atom_data", {}).items()
        },
    }


def build(state, representation, color, p, data, quality):
    p = dict(p)
    if color == "uniform":
        p["color"] = p.get("colorParams", {}).get("value", 0x808080)
    if representation in PARTICLE:
        names = {
            "particle-entity": "entity",
            "particle-compartment": "compartment",
            "particle-hierarchy": "hierarchy",
            "particle-index": "index",
            "particle-attribute": "attribute",
        }
        if color in names:
            p["colorBy"] = (
                p.get("colorParams", {}).get("attribute", "value")
                if color == "particle-attribute"
                else names[color]
            )
            p["scalarColor"] = color == "particle-attribute"
        elif color not in ("auto", "keep", "uniform"):
            raise ValueError(f"Theme {color} is not applicable to particles")
    if representation in VOLUME and color not in (
        "auto",
        "keep",
        "uniform",
        "volume-value",
        "volume-instance",
        "volume-segment",
    ):
        raise ValueError(f"Theme {color} is not applicable to a volume")
    if representation in STRUCTURE:
        if state is None or not state.atoms:
            raise ValueError(f"{representation} requires molecular atoms")
        return molecule.geometry(state, representation, color, p, data, quality)
    if representation in VOLUME:
        if "grid" not in data:
            raise ValueError(f"{representation} requires a local map or volume object")
        return volume.volume_geometry(data["grid"], representation, p, quality)
    if representation in PARTICLE:
        return shapes.particles(data, representation, p)
    if representation in SHAPE or representation in (
        "custom-label",
        "annotation-label",
    ):
        return (
            shapes.unitcell(state, p, data)
            if representation == "unitcell"
            else shapes.measurements(state, representation, p, data)
        )
    if representation in ("interactions", "cross-link-restraint", "clashes"):
        if state is None:
            raise ValueError("Contact geometry requires a structure")
        return analysis.geometry(state, representation, p, data)
    if representation in (
        "membrane-orientation",
        "assembly-symmetry",
        "confal-pyramids",
        "ntc-tube",
        "tunnel",
    ):
        return shapes.annotated(state, representation, p, data)
    if representation in ("orbital", "orbital-density"):
        return orbitals.geometry(data, representation, p)
    if representation == "mesh":
        return shapes.mesh_geometry(data, p)
    if representation == "g3d":
        return formats.g3d(data, p)
    if representation == "kinemage":
        return formats.kinemage(data, p)
    if representation == "mvs":
        return formats.mvs(data, p, quality, build)
    if representation == "pairwise-metric":
        g = Geometry()
        g.panels.append(panels.metric(data, state.atoms if state else ()))
        return g
    if representation not in PRESETS:
        raise ValueError(f"Unknown representation {representation}")
    if representation == "empty":
        return Geometry()
    if state is None or not state.atoms:
        raise ValueError(f"{representation} requires molecular atoms")
    if representation in ("auto", "auto-lod", "mesoscale"):
        count = len(state.atoms)
        rep = (
            "polymer-and-ligand"
            if count < 50_000
            else "polymer-cartoon"
            if count < 250_000
            else "coarse-surface"
        )
        return build(state, rep, color, p, data, quality)
    if (
        representation.startswith(("validation-", "quality-"))
        or representation == "partial-charges"
    ):
        mode = {
            "validation-geometry": "geometry-quality",
            "validation-density": "density-fit",
            "validation-rci": "random-coil-index",
            "quality-plddt": "plddt",
            "quality-qmean": "qmean",
            "partial-charges": "sb-ncbr-partial-charges",
        }[representation]
        return build(
            state,
            "polymer-and-ligand",
            mode if color == "auto" else color,
            p,
            data,
            quality,
        )
    if representation in (
        "atomic-detail",
        "illustrative",
        "coarse-surface",
        "polymer-cartoon",
    ):
        rep = {
            "atomic-detail": "ball-and-stick",
            "illustrative": "spacefill",
            "coarse-surface": "gaussian-surface",
            "polymer-cartoon": "cartoon",
        }[representation]
        return build(
            state,
            rep,
            "illustrative"
            if representation == "illustrative" and color == "auto"
            else color,
            p,
            data,
            quality,
        )
    result = Geometry()
    polymer = np.array([a.kind in ("protein", "nucleic") for a in state.atoms])
    carbohydrate = (
        np.array([a.resn in tables()["saccharides"] for a in state.atoms]) & ~polymer
    )
    water = np.array([a.resn in ("HOH", "WAT", "DOD") for a in state.atoms]) & ~polymer
    other = ~polymer & ~carbohydrate & ~water
    layers = [
        (polymer, "cartoon", 1),
        (other, "ball-and-stick", 1),
        (water, "ball-and-stick", 0.6),
    ]
    if representation == "protein-and-nucleic":
        layers = [(polymer, "cartoon", 1)]
    else:
        layers += [
            (carbohydrate, "carbohydrate", 1),
            (carbohydrate, "ball-and-stick", 0.3),
        ]
    for mask, rep, alpha in layers:
        if not mask.any():
            continue
        g = build(state.subset(mask), rep, color, p, subset_data(data, mask), quality)
        for piece in g.pieces:
            piece.mesh.opacity *= alpha
        result.extend(g)
    return result
