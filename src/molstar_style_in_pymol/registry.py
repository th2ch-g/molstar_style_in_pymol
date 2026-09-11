"""Public vocabulary and immutable reference revision; no PyMOL import."""

import json
from copy import deepcopy
from dataclasses import dataclass, field
from importlib.resources import files

REFERENCE_REVISION = "5b1b54ed03b03936041f514b33b8bb129b774d37"
STRUCTURE = (
    "cartoon",
    "backbone",
    "ball-and-stick",
    "blob-surface",
    "carbohydrate",
    "ellipsoid",
    "gaussian-surface",
    "gaussian-volume",
    "label",
    "line",
    "molecular-surface",
    "orientation",
    "plane",
    "point",
    "putty",
    "spacefill",
    "polyhedron",
)
VOLUME = ("direct-volume", "dot", "isosurface", "segment", "slice")
PARTICLE = (
    "particle-spacefill",
    "particle-orientation",
    "particle-fibers",
    "particle-target",
)
SHAPE = (
    "distance",
    "angle",
    "dihedral",
    "shape-label",
    "shape-orientation",
    "shape-plane",
    "unitcell",
)
EXTENSION = (
    "interactions",
    "cross-link-restraint",
    "membrane-orientation",
    "assembly-symmetry",
    "confal-pyramids",
    "ntc-tube",
    "clashes",
    "orbital",
    "orbital-density",
    "tunnel",
    "mesh",
    "kinemage",
    "g3d",
    "mvs",
    "pairwise-metric",
    "annotation-label",
    "custom-label",
)
PRESETS = (
    "default",
    "auto",
    "empty",
    "polymer-and-ligand",
    "protein-and-nucleic",
    "polymer-cartoon",
    "atomic-detail",
    "coarse-surface",
    "illustrative",
    "auto-lod",
    "mesoscale",
    "validation-geometry",
    "validation-density",
    "validation-rci",
    "quality-plddt",
    "quality-qmean",
    "partial-charges",
)
MATERIALS = {
    "matte": {"metalness": 0.0, "roughness": 1.0, "bumpiness": 0.0},
    "plastic": {"metalness": 0.0, "roughness": 0.2, "bumpiness": 0.0},
    "glossy": {"metalness": 0.0, "roughness": 0.6, "bumpiness": 0.0},
    "metallic": {"metalness": 1.0, "roughness": 0.6, "bumpiness": 0.0},
}
EFFECT_STYLES = (
    "outline",
    "occlusion",
    "shadow",
    "cel",
    "xray",
    "unlit",
    "bloom",
    "dof",
    "illumination",
    "background",
    "antialias",
)
OPERATIONS = ("list", "help", "refresh", "reset", "png", "ray")
STYLES = tuple(
    dict.fromkeys(
        (
            *PRESETS,
            *STRUCTURE,
            *VOLUME,
            *PARTICLE,
            *SHAPE,
            *EXTENSION,
            *MATERIALS,
            *EFFECT_STYLES,
        )
    )
)
QUALITIES = {
    "lowest": (2, 6),
    "lower": (3, 6),
    "low": (4, 8),
    "medium": (8, 12),
    "high": (12, 20),
    "higher": (16, 24),
    "highest": (24, 32),
    "auto": (8, 12),
    "custom": (8, 12),
}
ALIASES = {
    "cpk": "spacefill",
    "ballstick": "ball-and-stick",
    "ribbon": "cartoon",
    "surface": "molecular-surface",
    "sticks": "ball-and-stick",
    "tube": "backbone",
    "nucleic": "cartoon",
    "pae": "pairwise-metric",
    "segmentation": "segment",
    "tunnels": "tunnel",
    "snfg": "carbohydrate",
    "plddt": "quality-plddt",
    "qmean": "quality-qmean",
    "confal": "confal-pyramids",
}


def canonical(value):
    value = str(value).strip().lower().replace("_", "-")
    return ALIASES.get(value, value)


def reference():
    return json.loads(files(__package__).joinpath("reference.json").read_text())


@dataclass
class Profile:
    representation: str
    params: dict = field(default_factory=dict)
    material: dict = field(default_factory=lambda: deepcopy(MATERIALS["matte"]))
    effects: dict = field(default_factory=dict)
    edges: str = "none"
    edge_width: float = 0.04
    ignore_light: bool = False
    background: object = None


def resolve(style, representation, params):
    style = canonical(style)
    if style not in STYLES:
        raise ValueError(f"Unknown style {style!r}; use molstar_style list")
    rep = (
        canonical(representation) if representation not in (None, "", "auto") else style
    )
    p = deepcopy(params)
    material = p.pop("material", "matte")
    if isinstance(material, str):
        if material not in MATERIALS:
            raise ValueError(
                f"material must be one of {tuple(MATERIALS)} or a dictionary"
            )
        material = MATERIALS[material]
    material = {**MATERIALS["matte"], **material}
    if set(material) - set(MATERIALS["matte"]) or any(
        not 0 <= float(v) <= 1 for v in material.values()
    ):
        raise ValueError("material accepts metalness, roughness, bumpiness in [0, 1]")
    effects = deepcopy(p.pop("postprocessing", {}))
    if style in MATERIALS:
        material = deepcopy(MATERIALS[style])
        rep = "polymer-and-ligand" if rep == style else rep
    if style in EFFECT_STYLES:
        if style not in ("unlit", "cel", "xray", "background"):
            effects.setdefault(style, True)
        else:
            p[style] = True
        rep = "polymer-and-ligand" if rep == style else rep
    if style == "background":
        p.setdefault(
            "background", {"variant": {"name": "horizontalGradient", "params": {}}}
        )
    if p.get("illumination"):
        effects["illumination"] = p["illumination"]
    if rep == "default":
        rep = "polymer-and-ligand"
    if rep not in STYLES:
        raise ValueError(f"Unknown representation {rep!r}")
    # Mol* enables occlusion by default; it is scoped to the managed layer.
    effects.setdefault("occlusion", True)
    if style == "illustrative":
        p.setdefault("ignoreLight", True)
        effects.setdefault("outline", True)
    return Profile(
        rep,
        p,
        material,
        effects,
        "silhouette" if effects.get("outline") else "none",
        float(p.pop("edgeWidth", 0.04)),
        bool(p.get("ignoreLight", p.get("unlit", False))),
    )
