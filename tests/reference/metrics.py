"""Compare registered renders without fitting colors, masks, or image alignment."""

import argparse
import json
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.ndimage import binary_erosion


def compare(reference, actual):
    a = np.asarray(Image.open(reference).convert("RGB"), dtype=float)
    b = np.asarray(Image.open(actual).convert("RGB"), dtype=float)
    if a.shape != b.shape:
        raise ValueError("Images must have the same dimensions")
    foreground_a, foreground_b = a.min(axis=2) < 235, b.min(axis=2) < 235
    intersection = foreground_a & foreground_b
    union = foreground_a | foreground_b
    interior = binary_erosion(intersection, iterations=2)
    if not interior.any():
        raise ValueError("No shared foreground interior")
    return {
        "silhouette_iou": float(intersection.sum() / union.sum()),
        "interior_rgb_mae_255": float(abs(a - b)[interior].mean()),
        "interior_pixels": int(interior.sum()),
        "width": a.shape[1],
        "height": a.shape[0],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("actual", type=Path)
    args = parser.parse_args()
    print(json.dumps(compare(args.reference, args.actual), indent=2))
