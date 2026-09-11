"""Mol* alpha-orbital spherical Gaussian evaluation using local basis data."""

import numpy as np

from .data import finite
from .harmonics import FUNCTIONS
from .scene import Geometry, PickTarget
from .themes import rgb
from .volume import grid_box, isosurface


def grid(data, p, density=False):
    basis = data.get("basis", {}).get("atoms")
    orbitals = data.get("orbitals")
    if not basis or not orbitals:
        raise ValueError("Orbital computation requires local basis.atoms and orbitals")
    coords = finite([a["center"] for a in basis], (None, 3), "basis centers in Bohr")
    g = grid_box(
        coords,
        np.full(len(coords), float(p.get("boxExpand", 4.5))),
        float(p.get("gridSpacing", 0.4)),
        int(p.get("_budget", 512 * 1024**2)),
    )
    points = g.world(np.indices(g.values.shape).reshape(3, -1).T)
    selected = orbitals if density else [orbitals[int(p.get("index", 0))]]
    total = np.zeros(len(points), np.float32)
    order = p.get("sphericalOrder", data.get("sphericalOrder", "gaussian"))
    if order not in ("gaussian", "cca", "cca-reverse"):
        raise ValueError("sphericalOrder must be gaussian, cca, or cca-reverse")
    for orbital in selected:
        value = np.zeros(len(points))
        offset = 0
        alpha = finite(orbital.get("alpha"), (None,), "orbital alpha coefficients")
        for atom in basis:
            x, y, z = (points - np.asarray(atom["center"])).T
            r2 = x * x + y * y + z * z
            for shell in atom["shells"]:
                exponents = finite(shell["exponents"], (None,), "Gaussian exponents")
                if np.any(exponents <= 0):
                    raise ValueError("Gaussian exponents must be positive")
                for j, degree_index in enumerate(shell["angularMomentum"]):
                    if not isinstance(degree_index, int) or not 0 <= degree_index <= 4:
                        raise ValueError("Mol* spherical basis supports L=0..4")
                    coeff = finite(
                        shell["coefficients"][j],
                        (len(exponents),),
                        "contraction coefficients",
                    )
                    weights = alpha[offset : offset + 2 * degree_index + 1]
                    offset += 2 * degree_index + 1
                    if len(weights) != 2 * degree_index + 1:
                        raise ValueError("Orbital alpha length does not match basis")
                    if order != "gaussian":
                        indices = [degree_index] + [
                            v
                            for m in range(1, degree_index + 1)
                            for v in (
                                (degree_index + m, degree_index - m)
                                if order == "cca"
                                else (degree_index - m, degree_index + m)
                            )
                        ]
                        weights = weights[indices]
                    radial = sum(c * np.exp(-e * r2) for c, e in zip(coeff, exponents))
                    cutoff = float(p.get("cutoffThreshold", 0.001))
                    if cutoff > 0:
                        radial[r2 > -np.log(cutoff) / exponents.min()] = 0
                    value += radial * FUNCTIONS[degree_index](weights, x, y, z)
        if offset != len(alpha):
            raise ValueError("Unused orbital alpha coefficients")
        total += (
            float(orbital.get("occupancy", 0)) * value * value if density else value
        )
    g.values = total.reshape(g.values.shape)
    g.transform[:3, :] *= 0.529177210859
    g.label = "electron density" if density else f"orbital {p.get('index', 0)}"
    return g


def geometry(data, rep, p):
    g = data.get("grid") or grid(data, p, rep == "orbital-density")
    result = Geometry()
    target = (PickTarget(label=g.label),)
    threshold = float(
        p.get("isoValue", max(abs(g.values.min()), abs(g.values.max())) * 0.15)
    )
    if threshold <= 0:
        raise ValueError("Orbital isoValue must be positive")
    if g.values.max() > threshold:
        result.add(
            isosurface(g, threshold, rgb(p.get("positiveColor", 0x3366FF))), target
        )
    if rep == "orbital" and g.values.min() < -threshold:
        result.add(
            isosurface(g, -threshold, rgb(p.get("negativeColor", 0xFF3333))), target
        )
    if not result.pieces:
        raise ValueError("Orbital isovalue is outside the computed field")
    return result
