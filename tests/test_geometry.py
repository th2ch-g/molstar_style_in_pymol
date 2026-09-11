"""Coverage checks assert physical geometry, ownership, and reference formulas."""

import numpy as np
import pytest
from samples import data_for, molecule, scalar_grid

from molstar_style_in_pymol.analysis import accessible_area
from molstar_style_in_pymol.geometry import build
from molstar_style_in_pymol.registry import (
    EXTENSION,
    PARTICLE,
    PRESETS,
    SHAPE,
    STRUCTURE,
    VOLUME,
    reference,
)
from molstar_style_in_pymol.volume import density_slices, isosurface


@pytest.fixture(scope="module")
def state():
    return molecule()


@pytest.mark.parametrize(
    "rep",
    (
        *STRUCTURE,
        *VOLUME,
        *PARTICLE,
        *SHAPE,
        *(r for r in EXTENSION if r != "mvs"),
        *(r for r in PRESETS if r != "empty"),
    ),
)
def test_geometry(rep, state):
    g = build(state, rep, "auto", {"resolution": 1.2}, data_for(rep, state), "low")
    assert g.pieces or g.volumes or g.panels, rep
    for piece in g.pieces:
        mesh = piece.mesh
        assert len(mesh.faces) > 0
        assert np.isfinite(mesh.vertices).all() and np.isfinite(mesh.normals).all()
        assert np.isfinite(mesh.colors).all()
        assert mesh.faces.max() < len(mesh.vertices)
        assert mesh.owners.min() >= 0 and mesh.owners.max() < len(piece.atoms)
        assert np.linalg.norm(mesh.normals, axis=1).max() <= 1.001


@pytest.mark.parametrize(
    "rep,visual",
    [(r, v) for r, values in reference()["visuals"].items() for v in values],
)
def test_every_visual(rep, visual, state):
    g = build(
        state,
        rep,
        "auto",
        {"visuals": [visual], "resolution": 1.2},
        data_for(rep, state),
        "low",
    )
    assert g.pieces or g.volumes, (rep, visual)


def test_density_samples_are_transparent_planes():
    g = build(None, "direct-volume", "auto", {}, data_for("direct-volume"), "low")
    for direction in (None, [0, 0, 1], [1, 0.2, 0.3]):
        mesh = density_slices(g.volumes[0], "low", direction)
        assert len(mesh.faces) > 0
        assert np.any((mesh.alphas > 0) & (mesh.alphas < 1))


def test_isosurface_coordinate_transform():
    g = scalar_grid()
    m = isosurface(g, 0.5)
    np.testing.assert_allclose(
        m.vertices.mean(axis=0), g.world([[5, 6, 7]])[0], atol=0.1
    )
    np.testing.assert_allclose(np.linalg.norm(m.normals, axis=1), 1, atol=1e-5)


def test_isolated_atom_asa(state):
    value = accessible_area(state.atoms[:1], state.coords[:1])
    np.testing.assert_allclose(value, [4 * np.pi * (state.atoms[0].vdw + 1.4) ** 2])


def test_measurement_values():
    for rep, expected in (("distance", 3), ("angle", 90), ("dihedral", 116.565051177)):
        g = build(None, rep, "auto", {}, data_for(rep), "low")
        assert abs(abs(g.metadata["value"]) - expected) < 1e-6


def test_orbital_p_symmetry_and_density():
    from molstar_style_in_pymol.orbitals import grid

    data = data_for("orbital")
    g = grid(data, {"gridSpacing": 0.5})
    np.testing.assert_allclose(g.values, -g.values[:, :, ::-1], atol=1e-6)
    density = grid(data, {"gridSpacing": 0.5}, True)
    np.testing.assert_allclose(density.values, 2 * g.values**2, atol=1e-6)


def test_unknown_theme_and_missing_annotation_fail(state):
    with pytest.raises(ValueError, match="requires"):
        build(state, "quality-plddt", "auto", {}, {}, "low")
    with pytest.raises(ValueError, match="Unknown color"):
        build(state, "cartoon", "made-up", {}, {}, "low")
