"""Local geometric analysis; annotations remain distinct from predictions."""

import numpy as np
from scipy.spatial import cKDTree

from .mesh import merge, unit
from .primitives import dashed, sphere
from .scene import Geometry
from .themes import rgb

INTERACTION_COLORS = {
    "hydrogen-bond": 0x2B83BA,
    "weak-hydrogen-bond": 0x92C5DE,
    "hydrophobic": 0x808080,
    "ionic": 0xF6C343,
    "metal-coordination": 0x8C6BB1,
    "halogen-bond": 0xD7191C,
    "pi-stacking": 0x66C2A5,
    "cation-pi": 0xFDAE61,
    "clash": 0xFF0000,
}


def accessible_area(atoms, coords, params=None):
    p = params or {}
    count = int(p.get("numberOfSpherePoints", 92))
    probe = float(p.get("probeSize", 1.4))
    if count < 12 or count > 10000 or probe < 0:
        raise ValueError(
            "ASA requires 12..10000 sphere samples and nonnegative probeSize"
        )
    i = np.arange(count)
    z = 1 - 2 * (i + 0.5) / count
    a = i * np.pi * (3 - np.sqrt(5))
    radius = np.sqrt(1 - z * z)
    points = np.c_[np.cos(a) * radius, np.sin(a) * radius, z]
    radii = np.array([v.vdw for v in atoms]) + probe
    tree = cKDTree(coords)
    areas = np.zeros(len(atoms))
    for i, (center, r) in enumerate(zip(coords, radii)):
        neighbors = [
            j for j in tree.query_ball_point(center, r + radii.max()) if j != i
        ]
        exposed = np.ones(count, bool)
        samples = center + points * r
        for j in neighbors:
            exposed &= np.sum((samples - coords[j]) ** 2, axis=1) >= radii[j] ** 2
        areas[i] = 4 * np.pi * r * r * exposed.mean()
    return areas


def interactions(state, p):
    atoms, xyz = state.atoms, state.coords
    neighbors = [set() for _ in atoms]
    for i, j in state.bonds:
        neighbors[i].add(j)
        neighbors[j].add(i)
    pairs = cKDTree(xyz).query_pairs(
        float(p.get("maxDistance", 5.5)), output_type="ndarray"
    )
    records = []
    metals = {"ZN", "FE", "MG", "MN", "CA", "CU", "CO", "NI", "CD"}

    def charge(a):
        if a.formal_charge:
            return a.formal_charge
        if a.resn in ("LYS", "ARG") and a.element == "N":
            return 1
        if a.resn in ("ASP", "GLU") and a.name in ("OD1", "OD2", "OE1", "OE2"):
            return -1
        return 0

    def donor(i):
        a = atoms[i]
        h = [j for j in neighbors[i] if atoms[j].element in ("H", "D")]
        return h, a.element in ("N", "O", "S") and (
            bool(h)
            or (a.element == "N" and a.formal_charge >= 0)
            or a.resn in ("SER", "THR", "TYR", "HOH", "CYS")
        )

    for i, j in pairs:
        a, b = atoms[i], atoms[j]
        if (
            j in neighbors[i]
            or neighbors[i] & neighbors[j]
            or (a.model, a.chain, a.resi, a.segi) == (b.model, b.chain, b.resi, b.segi)
        ):
            continue
        if a.element in ("H", "D") or b.element in ("H", "D"):
            continue
        distance = np.linalg.norm(xyz[i] - xyz[j])
        kind = None
        if (
            a.element in metals
            and b.element in ("O", "N", "S")
            or b.element in metals
            and a.element in ("O", "N", "S")
        ) and distance < 3:
            kind = "metal-coordination"
        elif charge(a) * charge(b) < 0 and distance < 5:
            kind = "ionic"
        elif (
            a.element in ("N", "O", "S")
            and b.element in ("N", "O", "S")
            and distance < 3.5
        ):
            for d, acceptor in ((i, j), (j, i)):
                h, isdonor = donor(d)
                if atoms[acceptor].formal_charge > 0 or not isdonor:
                    continue
                if h and not any(
                    np.dot(unit(xyz[d] - xyz[k]), unit(xyz[acceptor] - xyz[k])) < -0.5
                    for k in h
                ):
                    continue
                kind = "hydrogen-bond"
                break
        elif (
            (a.element in ("CL", "BR", "I") and b.element in ("N", "O", "S"))
            or (b.element in ("CL", "BR", "I") and a.element in ("N", "O", "S"))
        ) and distance < 4:
            halogen, acceptor = (i, j) if a.element in ("CL", "BR", "I") else (j, i)
            if any(
                np.dot(unit(xyz[k] - xyz[halogen]), unit(xyz[acceptor] - xyz[halogen]))
                < -0.7
                for k in neighbors[halogen]
            ):
                kind = "halogen-bond"
        elif a.element == "C" and b.element == "C" and distance < 4:
            kind = "hydrophobic"
        elif (
            {a.element, b.element} <= {"C", "O", "N"}
            and distance < 3.8
            and p.get("weakHydrogenBonds", False)
        ):
            kind = "weak-hydrogen-bond"
        if kind:
            records.append(
                {"indices": [int(i), int(j)], "type": kind, "distance": float(distance)}
            )
    from .molecule import axes_of, residues

    rings = []
    for ids in residues(state):
        ring = [
            i
            for i in ids
            if atoms[i].resn in ("PHE", "TYR", "TRP", "HIS")
            and atoms[i].name not in ("N", "CA", "C", "O", "CB")
            and atoms[i].element in ("C", "N")
        ]
        if len(ring) >= 5:
            center, axes, _ = axes_of(xyz[ring])
            rings.append((ring[0], center, axes[2]))
    for r, (i, a, an) in enumerate(rings):
        for j, b, bn in rings[r + 1 :]:
            if np.linalg.norm(a - b) < 5.5 and (
                abs(np.dot(an, bn)) > 0.7 or abs(np.dot(an, bn)) < 0.3
            ):
                records.append(
                    {"indices": [i, j], "positions": [a, b], "type": "pi-stacking"}
                )
        for j, atom in enumerate(atoms):
            if (
                charge(atom) > 0
                and np.linalg.norm(xyz[j] - a) < 5
                and abs(np.dot(unit(xyz[j] - a), an)) > 0.7
            ):
                records.append(
                    {"indices": [i, j], "positions": [a, xyz[j]], "type": "cation-pi"}
                )
    return records


def geometry(state, representation, p, data):
    records = data.get(
        "interactions"
        if representation == "interactions"
        else "cross_links"
        if representation == "cross-link-restraint"
        else "clashes"
    )
    if records is None and representation == "interactions":
        records = interactions(state, p)
    if records is None and representation == "clashes" and p.get("compute", False):
        bonded = {tuple(sorted(b)) for b in state.bonds}
        records = [
            {"indices": [int(i), int(j)], "type": "clash"}
            for i, j in cKDTree(state.coords).query_pairs(4)
            if (i, j) not in bonded
            and np.linalg.norm(state.coords[i] - state.coords[j])
            < state.atoms[i].vdw + state.atoms[j].vdw - float(p.get("overlap", 0.4))
        ]
    if records is None:
        raise ValueError(
            f"{representation} requires local annotations; computed clashes need compute=true"
        )
    pieces = []
    for row in records:
        indices = row.get("indices")
        if (
            indices is None
            or len(indices) != 2
            or any(int(i) != i or not 0 <= i < len(state.atoms) for i in indices)
        ):
            raise ValueError(
                "Contact indices must be two zero-based selected atom positions"
            )
        i, j = indices
        a, b = np.asarray(row.get("positions", state.coords[[i, j]]), float)
        kind = row.get(
            "type", "clash" if representation == "clashes" else "hydrogen-bond"
        )
        if p.get("types") and kind not in p["types"]:
            continue
        color = rgb(row.get("color", INTERACTION_COLORS.get(kind, 0xFF8800)))
        if representation == "cross-link-restraint":
            distance = float(np.linalg.norm(a - b))
            lower = float(row.get("lower", 0))
            upper = float(row.get("upper", np.inf))
            color = rgb(0xFF3333 if distance < lower or distance > upper else 0x33AA55)
        pieces.append(dashed(a, b, float(p.get("sizeFactor", 0.08)), color, i))
        if representation == "clashes":
            pieces.append(
                sphere((a + b) / 2, float(p.get("sizeFactor", 0.25)), color, i, 8)
            )
    result = Geometry().add(merge(pieces), state.atoms)
    result.metadata["contacts"] = records
    return result
