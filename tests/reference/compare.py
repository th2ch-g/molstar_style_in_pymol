"""Compare Python curves with values executed by the pinned Mol* checkout."""

import argparse
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from molstar_style_in_pymol.molecule import polymer, residues
from molstar_style_in_pymol.polymer_trace import (
    controls_for,
    helix_centers,
    interpolate,
)
from molstar_style_in_pymol.source import read


def compare(cmd, structure, reference):
    cmd.load(str(structure), "sample")
    snapshots, _ = read(cmd, "sample and polymer")
    state = next(iter(snapshots.values()))[0]
    records = {(s["chain"], s["resi"]): s for s in reference["segments"]}
    native_mismatch = 0
    groups = residues(state)
    chains = {}
    for group in groups:
        names = {
            state.atoms[i].name: i for i in group if state.atoms[i].alt in ("", "A")
        }
        atom = state.atoms[group[0]]
        key = (atom.chain, atom.resi)
        if key not in records:
            continue
        anchor = records[key].get("trace", "CA")
        if anchor not in names:
            continue
        atom = state.atoms[names[anchor]]
        native_mismatch += atom.ss != records[key]["ss"]
        point = state.coords[names[anchor]]
        pair = (
            (("C3'", "C1'") if atom.resn.startswith("D") else ("C4'", "C3'"))
            if atom.kind == "nucleic"
            else ("C", "O")
        )
        direction = state.coords[names[pair[1]]] - state.coords[names[pair[0]]]
        chains.setdefault(atom.chain, []).append((point, direction, records[key]))
    differences = {"control": [], "curve": [], "normal": [], "binormal": []}
    for rows in chains.values():
        chunks = []
        for row in rows:
            if chunks and np.linalg.norm(row[0] - chunks[-1][-1][0]) < (
                10 if row[2].get("kind") == "nucleic" else 4.5
            ):
                chunks[-1].append(row)
            else:
                chunks.append([row])
        for chunk in chunks:
            points = np.array([r[0] for r in chunk])
            directions = np.array([r[1] for r in chunk])
            secondary = [r[2]["ss"] for r in chunk]
            tubular = (
                reference.get("options", {})
                .get("params", {})
                .get("tubularHelices", False)
            )
            if tubular:
                centers = helix_centers(points)
                mask = np.array(secondary) == "H"
                points[mask] = centers[mask]
                directions[mask] = [1, 0, 0]
            for i, (_, _, record) in enumerate(chunk):
                control = controls_for(points, directions, secondary, i, tubular)
                for k in ("p0", "p1", "p2", "p3", "p4", "d12", "d23"):
                    differences["control"].extend(
                        np.ravel(np.asarray(control[k]) - record["controls"][k])
                    )
                curve, normals, binormals = interpolate(
                    control,
                    0.9 if secondary[i] == "H" and not tubular else 0.5,
                    0.3 if record.get("kind") == "nucleic" else 0.5,
                )
                for kind, array, ref in (
                    ("curve", curve, "curve"),
                    ("normal", normals, "normals"),
                    ("binormal", binormals, "binormals"),
                ):
                    differences[kind].extend(array.ravel() - record[ref])
    report = {
        k: {
            "max": float(np.max(np.abs(v))),
            "rms": float(np.sqrt(np.mean(np.square(v)))),
        }
        for k, v in differences.items()
    }
    report["native_secondary_structure_mismatches"] = native_mismatch
    report["residues"] = len(records)
    state = state.subset(np.array([(a.chain, a.resi) in records for a in state.atoms]))
    state.atoms = tuple(
        replace(a, ss=records[(a.chain, a.resi)]["ss"]) for a in state.atoms
    )
    params = {
        "linearSegments": 8,
        "radialSegments": 16,
        "visuals": ["polymer-trace"],
        **reference.get("options", {}).get("params", {}),
    }
    mesh = polymer(
        state,
        "cartoon",
        np.ones((len(state.atoms), 3)),
        np.ones(len(state.atoms)),
        params,
        "medium",
    )
    vertices = np.concatenate(
        [
            np.asarray(m["vertices"]).reshape(-1, 3)[: m["vertexCount"]]
            for m in reference["meshes"]
        ]
    )
    # Tube cap centers can exist at turn boundaries that PyMOL's H/S/coil
    # classification merges. Compare surface rings independently of these
    # internal disks and their triangulation-only center vertices.
    centers = cKDTree(
        np.concatenate(
            [np.asarray(s["curve"]).reshape(-1, 3) for s in reference["segments"]]
        )
    )
    reference_skin = vertices[centers.query(vertices)[0] > 1e-4]
    python_skin = mesh.vertices[centers.query(mesh.vertices)[0] > 1e-4]
    distances = np.r_[
        cKDTree(reference_skin).query(python_skin)[0],
        cKDTree(python_skin).query(reference_skin)[0],
    ]
    report["mesh_surface_vertex_distance_angstrom"] = {
        "max": float(distances.max()),
        "rms": float(np.sqrt(np.mean(distances**2))),
        "python_vertices": len(mesh.vertices),
        "reference_vertices": len(vertices),
    }
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--tolerance", type=float, default=1e-4)
    args = parser.parse_args()
    import pymol2

    with pymol2.PyMOL() as instance:
        report = compare(
            instance.cmd, args.structure, json.loads(args.reference.read_text())
        )
        print(json.dumps(report, indent=2))
        for key in (
            "control",
            "curve",
            "normal",
            "binormal",
            "mesh_surface_vertex_distance_angstrom",
        ):
            if report[key]["max"] > args.tolerance:
                raise SystemExit(f"{key} exceeds tolerance {args.tolerance}")
