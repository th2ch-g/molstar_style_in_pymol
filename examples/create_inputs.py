"""Create local volume, particle, annotation, and orbital inputs for the guide."""

import argparse
import json
from pathlib import Path

import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", type=Path, default=Path(".cache/examples"))
args = parser.parse_args()
args.output.mkdir(parents=True, exist_ok=True)
xyz = np.indices((32, 32, 32), dtype=float)
r2 = sum((v - 15.5) ** 2 for v in xyz)
values = np.exp(-r2 / 35)
transform = np.eye(4)
transform[:3, :3] *= 0.4
transform[:3, 3] = -6.2
np.savez_compressed(
    args.output / "density.npz", values=values.astype("f4"), transform=transform
)
np.savez_compressed(
    args.output / "segments.npz",
    values=np.where(values > 0.3, 2, np.where(values > 0.08, 1, 0)).astype("f4"),
    transform=transform,
)
inputs = {
    "particles": {
        "particles": [
            {
                "position": [0, 0, 0],
                "radius": 0.5,
                "entity": "A",
                "quaternion": [0, 0, 0, 1],
                "points": [[0, 0, 0], [2, 1, 0], [3, 2, 2]],
                "target": [4, 3, 2],
            }
        ]
    },
    "tunnel": {
        "positions": [[0, 0, 0], [2, 1, 0], [3, 2, 2]],
        "radii": [0.8, 0.4, 0.6],
    },
    "membrane": {
        "membrane": {
            "center": [0, 0, 0],
            "normal": [0, 0, 1],
            "thickness": 3,
            "radius": 5,
        }
    },
    "pae": {
        "matrix": [[0, 2, 4], [1, 0, 6], [3, 5, 0]],
        "title": "Synthetic pairwise metric",
    },
    "orbital": {
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
    },
    "density_params": {
        "transferFunction": [
            [0, "blue", 0],
            [0.1, "blue", 0],
            [0.4, "cyan", 0.15],
            [1, "yellow", 0.8],
        ]
    },
    "cartoon_params": {
        "tubularHelices": False,
        "material": "plastic",
        "postprocessing": {"occlusion": True, "outline": True},
    },
}
for name, data in inputs.items():
    (args.output / f"{name}.json").write_text(json.dumps(data, indent=2) + "\n")
print(f"Wrote local examples to {args.output}")
