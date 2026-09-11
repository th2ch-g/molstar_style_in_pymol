"""Mol* color and size themes with explicit annotation requirements."""

from functools import lru_cache

import numpy as np

from .data import atom_fields, finite
from .registry import reference

ANNOTATION_THEMES = (
    "accessible-surface-area",
    "interaction-type",
    "cross-link",
    "assembly-symmetry-cluster",
    "confal-pyramids",
    "ntc-tube",
    "geometry-quality",
    "density-fit",
    "random-coil-index",
    "pdbe-structure-quality-report",
    "plddt",
    "qmean",
    "sb-ncbr-partial-charges",
    "mvs-annotation",
    "mvs-multilayer",
    "g3d-haplotype",
    "g3d-chromosome",
    "g3d-region",
)
ALIASES = {
    "chain": "chain-id",
    "element": "element-symbol",
    "ss": "secondary-structure",
    "rainbow": "sequence-id",
    "bfactor": "uncertainty",
    "molstar": "auto",
}


@lru_cache(maxsize=1)
def tables():
    return reference()


def rgb(value):
    if isinstance(value, str):
        from PIL import ImageColor

        value = ImageColor.getrgb(value)
    if isinstance(value, (int, np.integer)):
        return (
            np.array(
                [(int(value) >> 16) & 255, (int(value) >> 8) & 255, int(value) & 255],
                float,
            )
            / 255
        )
    values = np.array(value, float)
    if values.shape != (3,) or not np.isfinite(values).all():
        raise ValueError("Color must be an RGB triple, integer, or color name")
    return np.clip(values / (255 if values.max(initial=0) > 1 else 1), 0, 1)


def palette(value="many-distinct"):
    if isinstance(value, dict):
        value = value.get("colors", value.get("name", "many-distinct"))
    if isinstance(value, str):
        if value not in tables()["palettes"]:
            raise ValueError(f"Unknown palette: {value}")
        value = tables()["palettes"][value]
    return np.array([rgb(v) for v in value])


def categorical(values, colors="many-distinct"):
    colors = palette(colors)
    ids = {}
    return np.array(
        [colors[ids.setdefault(str(v), len(ids)) % len(colors)] for v in values], float
    ).reshape(-1, 3)


def scale(values, colors="red-blue", domain=None, missing=0xCCCCCC):
    values = np.array([np.nan if v is None else float(v) for v in values])
    valid = np.isfinite(values)
    if not valid.any():
        return np.tile(rgb(missing), (len(values), 1))
    lo, hi = (
        domain if domain is not None else (values[valid].min(), values[valid].max())
    )
    colors = palette(colors)
    if len(colors) == 1:
        out = np.tile(colors[0], (len(values), 1))
        out[~valid] = rgb(missing)
        return out
    t = np.clip((values - lo) / (hi - lo if hi != lo else 1), 0, 1) * (len(colors) - 1)
    t[~valid] = 0
    index = np.minimum(t.astype(int), len(colors) - 2)
    fraction = t - index
    out = (
        colors[index] * (1 - fraction[:, None]) + colors[index + 1] * fraction[:, None]
    )
    out[~valid] = rgb(missing)
    return out


def adjust(colors, saturation=0, lightness=0):
    if not saturation and not lightness:
        return colors
    from skimage.color import lab2rgb, rgb2lab

    lab = rgb2lab(np.asarray(colors).reshape(-1, 1, 3)).reshape(-1, 3)
    chroma = np.linalg.norm(lab[:, 1:], axis=1)
    ratio = np.divide(
        np.maximum(0, chroma + 18 * saturation),
        chroma,
        out=np.ones_like(chroma),
        where=chroma > 1e-4,
    )
    lab[:, 1:] *= ratio[:, None]
    lab[:, 0] += 18 * lightness
    return np.clip(lab2rgb(lab.reshape(-1, 1, 3)).reshape(-1, 3), 0, 1)


def required(fields, key, count):
    if key not in fields:
        raise ValueError(
            f"Theme requires atom_data.{key} or matching residue annotations"
        )
    if len(fields[key]) != count:
        raise ValueError(f"Theme field {key!r} does not match the selected atoms")
    return fields[key]


def colors(
    atoms, coords, mode="auto", representation="cartoon", params=None, data=None
):
    params, data = params or {}, data or {}
    mode = ALIASES.get(mode, mode)
    fields = atom_fields(atoms, data)
    n = len(atoms)
    if mode == "auto":
        mode = (
            "element-symbol"
            if representation
            in (
                "ball-and-stick",
                "line",
                "point",
                "ellipsoid",
                "spacefill",
                "polyhedron",
                "plane",
            )
            else "chain-id"
        )
    if mode == "keep":
        return np.array([a.color for a in atoms]).reshape(-1, 3)
    if mode == "uniform":
        return np.tile(rgb(params.get("value", 0x808080)), (n, 1))
    if mode in ("cartoon", "illustrative"):
        style = params.get("style", {"name": "chain-id"})
        base = colors(
            atoms,
            coords,
            style.get("name", "chain-id"),
            representation,
            style.get("params", {}),
            data,
        )
        if mode == "illustrative":
            base = adjust(base, lightness=float(params.get("lightness", 0.8)))
            for i, a in enumerate(atoms):
                if a.element.upper() not in ("C", "H"):
                    base[i] = rgb(0xEEEEEE)
        return base
    if mode == "element-symbol":
        table = tables()["color_maps"]["ElementSymbolColors"]
        custom = params.get("colors", {})
        if isinstance(custom, dict) and custom.get("name") == "custom":
            table = {**table, **custom.get("params", {})}
        base = np.array([rgb(table.get(a.element.upper(), 0xFFFFFF)) for a in atoms])
        carbon = params.get("carbonColor", {"name": "chain-id"})
        if isinstance(carbon, str):
            carbon = {"name": carbon}
        if carbon.get("name") != "element-symbol":
            replacement = colors(
                atoms,
                coords,
                carbon.get("name", "chain-id"),
                representation,
                carbon.get("params", {}),
                data,
            )
            for i, a in enumerate(atoms):
                if a.element.upper() == "C":
                    base[i] = replacement[i]
        return adjust(
            base,
            float(params.get("saturation", 0)),
            float(params.get("lightness", 0.2)),
        )
    intrinsic = {
        "atom-id": [a.name for a in atoms],
        "chain-id": [a.chain for a in atoms],
        "model-index": [a.model for a in atoms],
        "structure-index": [a.model for a in atoms],
        "unit-index": [(a.model, a.segi, a.chain) for a in atoms],
        "polymer-id": [(a.model, a.chain, a.kind) for a in atoms],
        "polymer-index": [(a.model, a.chain, a.kind) for a in atoms],
        "element-index": list(range(n)),
        "trajectory-index": [data.get("state", 1)] * n,
    }
    category = {
        "entity-id": "entity_id",
        "entity-source": "entity_source",
        "operator-hkl": "operator_hkl",
        "operator-name": "operator_name",
        "shape-group": "group",
        "volume-instance": "instance",
        "volume-segment": "segment",
        "particle-compartment": "compartment",
        "particle-entity": "entity",
        "particle-hierarchy": "hierarchy",
        "particle-index": "particle_index",
    }
    if mode in intrinsic or mode in category:
        values = (
            intrinsic[mode]
            if mode in intrinsic
            else required(fields, category[mode], n)
        )
        value = (
            params.get("palette", {})
            .get("params", {})
            .get("list", params.get("palette", {}).get("colors", "many-distinct"))
            if isinstance(params.get("palette", {}), dict)
            else params["palette"]
        )
        return categorical(values, value)
    if mode in (
        "secondary-structure",
        "molecule-type",
        "residue-name",
        "residue-charge",
    ):
        names = {
            "secondary-structure": "SecondaryStructureColors",
            "molecule-type": "MoleculeTypeColors",
            "residue-name": "ResidueNameColors",
            "residue-charge": "ChargedResidueColors",
        }
        table = tables()["color_maps"].get(names[mode], {})
        keys = []
        for a in atoms:
            molecule = (
                "protein"
                if a.kind == "protein"
                else ("dna" if a.resn.startswith("D") else "rna")
                if a.kind == "nucleic"
                else "water"
                if a.resn in ("HOH", "WAT", "DOD")
                else "saccharide"
                if a.resn
                in data.get(
                    "saccharides",
                    ("NAG", "MAN", "BMA", "GLC", "BGC", "GAL", "FUC", "SIA"),
                )
                else "ion"
                if len(a.resn) <= 2
                else "other"
            )
            if mode == "secondary-structure":
                keys.append(
                    {"H": "alphaHelix", "S": "betaStrand"}.get(
                        a.ss,
                        "coil"
                        if molecule == "protein"
                        else "carbohydrate"
                        if molecule == "saccharide"
                        else molecule,
                    )
                )
            elif mode == "molecule-type":
                keys.append(molecule)
            else:
                keys.append(a.resn)
        base = np.array([rgb(table.get(k, 0x808080)) for k in keys])
        return adjust(
            base,
            float(params.get("saturation", -1 if mode == "secondary-structure" else 0)),
            float(params.get("lightness", 0)),
        )
    if mode == "carbohydrate-symbol":
        from .carbohydrate import sugar_info

        return np.array([rgb(sugar_info(a.resn)[1]) for a in atoms])
    if mode == "sequence-id":
        residues = {}
        values = [
            residues.setdefault((a.model, a.chain, a.resi), len(residues))
            for a in atoms
        ]
        return scale(values, params.get("list", "rainbow"))
    if mode == "hydrophobicity":
        index = {"DGwif": 0, "DGwoct": 1, "Oct-IF": 2}[params.get("scale", "DGwif")]
        table = tables()["hydrophobicity"]
        domain = [
            max(v[index] for v in table.values()),
            min(v[index] for v in table.values()),
        ]
        return scale(
            [table.get(a.resn, [0, 0, 0])[index] for a in atoms],
            params.get("list", "red-yellow-green"),
            domain,
        )
    scalar = {
        "occupancy": [a.occupancy for a in atoms],
        "uncertainty": [a.bfactor for a in atoms],
        "formal-charge": [a.formal_charge for a in atoms],
        "partial-charge": [a.partial_charge for a in atoms],
    }
    if mode in scalar:
        default_domain = (
            [0, 1]
            if mode == "occupancy"
            else [-3, 3]
            if mode == "formal-charge"
            else [-1, 1]
            if mode == "partial-charge"
            else [0, 100]
        )
        return scale(
            fields.get(mode.replace("-", "_"), scalar[mode]),
            params.get("list", "red-blue"),
            params.get("domain", default_domain),
        )
    if mode == "external-volume":
        grid = data.get("grid")
        if grid is None:
            raise ValueError("external-volume coloring requires a local grid")
        return scale(
            grid.sample(coords), params.get("list", "red-blue"), params.get("domain")
        )
    if mode == "external-structure":
        from scipy.spatial import cKDTree

        positions = finite(
            data.get("external_coordinates"), (None, 3), "external_coordinates"
        )
        values = np.array([rgb(v) for v in data.get("external_colors", [])])
        if len(values) != len(positions):
            raise ValueError(
                "external_coordinates and external_colors must have equal length"
            )
        distance, index = cKDTree(positions).query(coords)
        result = values[index]
        result[distance > float(params.get("maxDistance", np.inf))] = rgb(
            params.get("defaultColor", 0xCCCCCC)
        )
        return result
    if mode in ("particle-attribute", "volume-value"):
        values = required(fields, params.get("attribute", "value"), n)
        return scale(values, params.get("list", "viridis"), params.get("domain"))
    if mode == "accessible-surface-area":
        from .analysis import accessible_area

        values = accessible_area(atoms, coords, params)
        return scale(
            values,
            params.get(
                "list",
                "blue-white-red"
                if "blue-white-red" in tables()["palettes"]
                else "red-blue",
            ),
            params.get("domain"),
        )
    if mode in ANNOTATION_THEMES:
        key = mode.replace("-", "_")
        if mode == "mvs-multilayer":
            base = colors(atoms, coords, "uniform", params={"value": 0xFFFFFF})
            for layer in data.get("layers", []):
                mask = np.array(required(fields, layer["field"], n), bool)
                base[mask] = rgb(layer["color"])
            return base
        values = required(fields, key, n)
        if mode == "plddt":
            result = np.tile(rgb(0xAAAAAA), (n, 1))
            for i, v in enumerate(values):
                if v is not None:
                    result[i] = rgb(
                        0x0053D6
                        if float(v) > 90
                        else 0x65CBF3
                        if float(v) > 70
                        else 0xFFDB13
                        if float(v) > 50
                        else 0xFF7D45
                    )
            return result
        if mode in (
            "assembly-symmetry-cluster",
            "confal-pyramids",
            "ntc-tube",
            "interaction-type",
            "g3d-haplotype",
            "g3d-chromosome",
            "g3d-region",
            "mvs-annotation",
        ):
            return categorical(values, params.get("list", "many-distinct"))
        return scale(
            values, params.get("list", "red-yellow-green"), params.get("domain")
        )
    raise ValueError(f"Unknown color theme {mode!r}")


def sizes(atoms, mode="physical", params=None, data=None):
    params, data = params or {}, data or {}
    fields = atom_fields(atoms, data)
    n = len(atoms)
    if mode == "physical":
        result = np.array([a.vdw for a in atoms])
    elif mode == "uniform":
        result = np.full(n, float(params.get("value", 1)))
    elif mode == "uncertainty":
        result = float(params.get("baseSize", 0.2)) + np.array(
            [a.bfactor for a in atoms]
        ) * float(params.get("bfactorFactor", 0.1))
        if "rmsf" in fields:
            result = float(params.get("baseSize", 0.2)) + np.array(
                fields["rmsf"]
            ) * float(params.get("rmsfFactor", 0.05))
    elif mode in ("shape-group", "particle-size", "volume-value"):
        key = {
            "shape-group": "size",
            "particle-size": "radius",
            "volume-value": "value",
        }[mode]
        result = np.array(required(fields, key, n), float)
    else:
        raise ValueError(f"Unknown size theme {mode!r}")
    if not np.isfinite(result).all() or (result < 0).any():
        raise ValueError("Size theme produced invalid radii")
    return result
