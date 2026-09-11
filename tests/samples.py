"""Deterministic synthetic scientific inputs; no downloads or user files."""

import numpy as np
from chempy import Atom as ModelAtom
from chempy import Bond
from chempy.models import Indexed

from molstar_style_in_pymol.data import Grid
from molstar_style_in_pymol.source import Atom, State


def molecule():
    atoms = []
    coords = []
    bonds = []
    model = Indexed()

    def atom(
        name, pos, element="C", resn="ALA", resi="1", chain="A", kind="protein", ss="H"
    ):
        i = len(atoms)
        atoms.append(
            Atom(
                "sample",
                i + 1,
                name,
                resi,
                resn,
                chain,
                "",
                element,
                ss,
                kind,
                (0.2, 0.6, 0.9),
                1.52 if element == "O" else 1.7,
                bfactor=10 + i,
                aniso=(0.15, 0.22, 0.35, 0.02, 0.01, 0.02),
            )
        )
        coords.append(pos)
        a = ModelAtom()
        a.index = i + 1
        a.name = name
        a.resi = resi
        a.resn = resn
        a.chain = chain
        a.symbol = element
        a.coord = list(pos)
        a.ss = ss
        a.b = 10 + i
        a.vdw = atoms[-1].vdw
        a.hetatm = kind == "other"
        model.atom.append(a)
        return i

    def bond(i, j, order=1):
        b = Bond()
        b.index = [i, j]
        b.order = order
        model.bond.append(b)
        bonds.append((i, j))

    last = None
    for i in range(8):
        angle = i * np.deg2rad(100)
        center = np.array([2.3 * np.cos(angle), 2.3 * np.sin(angle), i * 1.5])
        if i >= 6:
            center += [7, 0, 0]
        ca = atom("CA", center, resi=str(i + 1), ss="H" if i < 4 else "S")
        carbon = atom(
            "C", center + [0.5, 0.4, 0.8], resi=str(i + 1), ss="H" if i < 4 else "S"
        )
        oxygen = atom(
            "O",
            center + [0.3, 1.3, 0.7],
            "O",
            resi=str(i + 1),
            ss="H" if i < 4 else "S",
        )
        bond(ca, carbon)
        bond(carbon, oxygen, 2)
        if last is not None and i != 6:
            bond(last, ca)
        last = carbon
    for k, center in enumerate((np.array([0, 12, 0]), np.array([4, 12, 0]))):
        group = []
        for i in range(6):
            t = i * np.pi / 3
            group.append(
                atom(
                    f"C{i + 1}",
                    center + [np.cos(t), np.sin(t), 0.15 * (-1) ** i],
                    resn="NAG" if k == 0 else "MAN",
                    resi=str(k + 1),
                    chain="C" if k == 0 else "D",
                    kind="other",
                    ss="",
                )
            )
        for i in range(6):
            bond(group[i], group[(i + 1) % 6])
        if k:
            bond(group[0] - 6, group[0])
    group = []
    for i, name in enumerate(("N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4")):
        t = i * 2 * np.pi / 9
        group.append(
            atom(
                name,
                [12 + 1.5 * np.cos(t), 1.5 * np.sin(t), 0],
                name[0],
                resn="DA",
                chain="B",
                kind="nucleic",
                ss="",
            )
        )
    for i in range(len(group)):
        bond(group[i], group[(i + 1) % len(group)])
    anchor = atom("C4'", [12, -3, 0], resn="DA", chain="B", kind="nucleic", ss="")
    bond(anchor, group[0])
    atom("P", [12, -4, 0], "P", resn="DA", chain="B", kind="nucleic", ss="")
    metal = atom("ZN", [0, -8, 0], "ZN", resn="ZN", chain="Z", kind="other", ss="")
    for k, v in enumerate(((1, 1, 1), (-1, -1, 1), (1, -1, -1), (-1, 1, -1))):
        j = atom(
            "O",
            np.array([0, -8, 0]) + np.array(v) * 1.15,
            "O",
            resn="HOH",
            resi=str(k + 1),
            chain="W",
            kind="other",
            ss="",
        )
        bond(metal, j)
    return State(
        tuple(atoms),
        np.asarray(coords, np.float32),
        tuple(bonds),
        model,
        np.full(len(atoms), 16),
        {},
    )


def scalar_grid(segment=False):
    xyz = np.indices((12, 13, 14), dtype=float)
    values = np.exp(-sum((xyz[i] - 5 - i) ** 2 for i in range(3)) / 8)
    if segment:
        values = np.where(values > 0.4, 2, np.where(values > 0.1, 1, 0))
    transform = np.eye(4)
    transform[:3, :3] *= 0.7
    transform[:3, 3] = [-4, -5, -6]
    return Grid(values, transform, "synthetic density")


def data_for(rep, state=None):
    if rep in ("direct-volume", "dot", "isosurface", "slice", "segment"):
        return {"grid": scalar_grid(rep == "segment")}
    if rep.startswith("particle-"):
        return {
            "particles": [
                {
                    "position": [0, 0, 0],
                    "radius": 0.5,
                    "target": [4, 3, 2],
                    "points": [[0, 0, 0], [2, 1, 0], [3, 2, 1]],
                    "quaternion": [0, 0, 0, 1],
                }
            ]
        }
    if rep in (
        "distance",
        "angle",
        "dihedral",
        "shape-plane",
        "shape-orientation",
        "shape-label",
        "annotation-label",
        "custom-label",
    ):
        n = {"distance": 2, "angle": 3, "dihedral": 4}.get(rep, 4)
        return {
            "positions": [[0, 0, 0], [3, 0, 0], [3, 2, 0], [4, 2, 2]][:n],
            "text": "Metric",
        }
    if rep == "unitcell":
        return {"cell": [10, 12, 14, 90, 100, 110]}
    if rep in ("interactions", "cross-link-restraint", "clashes"):
        return {
            {
                "interactions": "interactions",
                "cross-link-restraint": "cross_links",
                "clashes": "clashes",
            }[rep]: [
                {"indices": [0, 4], "type": "hydrogen-bond", "lower": 1, "upper": 3}
            ]
        }
    if rep == "membrane-orientation":
        return {
            "membrane": {
                "center": [0, 0, 0],
                "normal": [0, 0, 1],
                "thickness": 4,
                "radius": 7,
            }
        }
    if rep == "assembly-symmetry":
        return {
            "symmetry": {"axes": [{"start": [0, 0, -5], "end": [0, 0, 5], "order": 3}]}
        }
    if rep in ("confal-pyramids", "ntc-tube"):
        return {
            "steps": [
                {
                    "class": "AA00",
                    "score": 80,
                    "positions": [
                        [-2, 0, 0],
                        [0, 0, 0],
                        [0, 1, 1],
                        [0, -1, 1],
                        [2, 0, 0],
                    ],
                }
            ]
        }
    if rep == "tunnel":
        return {"positions": [[0, 0, 0], [2, 0, 1], [4, 1, 2]], "radii": [0.8, 0.5, 1]}
    if rep == "mesh":
        return {
            "vertices": [[0, 0, 0], [3, 0, 0], [0, 3, 0], [0, 0, 3]],
            "faces": [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]],
        }
    if rep in ("orbital", "orbital-density"):
        return {
            "basis": {
                "atoms": [
                    {
                        "center": [0, 0, 0],
                        "shells": [
                            {
                                "exponents": [0.5],
                                "angularMomentum": [1],
                                "coefficients": [[1]],
                            }
                        ],
                    }
                ]
            },
            "orbitals": [{"alpha": [1, 0, 0], "occupancy": 2, "energy": -0.5}],
        }
    if rep == "kinemage":
        return {
            "text": "@kinemage 1\n@vectorlist {trace} color=red\n{a} P 0 0 0\n{b} 3 2 1\n{c} 4 3 0\n@balllist {points} color=blue\n{d} 2 0 1"
        }
    if rep == "g3d":
        return {
            "haplotypes": {
                "maternal": {
                    "chr1": {
                        "start": [0, 1000, 2000],
                        "x": [0, 3, 4],
                        "y": [0, 2, 1],
                        "z": [0, 1, 3],
                    }
                }
            }
        }
    if rep == "pairwise-metric":
        n = (
            len({(a.model, a.segi, a.chain, a.resi) for a in state.atoms})
            if state
            else 4
        )
        return {
            "predicted_aligned_error": (
                np.abs(np.arange(n)[:, None] - np.arange(n)[None, :]) + 0.5
            ).tolist()
        }
    if rep.startswith(("validation-", "quality-")) or rep == "partial-charges":
        n = len(state.atoms)
        return {
            "atom_data": {
                k: np.linspace(0, 100, n)
                for k in (
                    "geometry_quality",
                    "density_fit",
                    "random_coil_index",
                    "plddt",
                    "qmean",
                    "sb_ncbr_partial_charges",
                )
            }
        }
    return {}
