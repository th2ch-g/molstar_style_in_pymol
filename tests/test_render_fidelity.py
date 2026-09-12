"""Regressions for reference ribbon geometry, physical materials and export."""

from types import SimpleNamespace

import numpy as np
import pytest

from molstar_style_in_pymol.lighting import settings, shade
from molstar_style_in_pymol.polymer_trace import controls_for, interpolate, section_mesh
from molstar_style_in_pymol.postprocessing import ssao_samples
from molstar_style_in_pymol.registry import resolve


def profile(material="matte", **lighting):
    return resolve(material, "cartoon", {"lighting": lighting})


def test_matte_uses_dielectric_fresnel_instead_of_white_phong_highlight():
    p = profile(
        ambientIntensity=0, light=[{"inclination": 180, "azimuth": 0, "intensity": 1}]
    )
    result = shade([[0.5, 0.5, 0.5]], [[0, 0, 1]], [[0, 0, 1]], p)
    # At normal incidence roughness=1: D=1/pi, V=1/4, diffuse=0.5.
    fresnel = 0.04 + 0.96 * 2 ** (-5.55473 - 6.98316)
    np.testing.assert_allclose(result, 0.5 + fresnel / 4, atol=1e-7)


def test_metal_has_no_diffuse_component():
    p = profile("metallic", ambientIntensity=0, light=[])
    np.testing.assert_allclose(
        shade([[0.2, 0.6, 0.4]], [[0, 0, 1]], [[0, 0, 1]], p), 0.01
    )


def test_reference_default_light_and_ambient():
    directions, colors, ambient, exposure = settings(profile())
    np.testing.assert_allclose(
        directions, [[-0.3830222216, 0.3213938048, 0.8660254038]], atol=1e-9
    )
    np.testing.assert_allclose(colors, 0.6)
    np.testing.assert_allclose(ambient, 0.4)
    assert exposure == 1


@pytest.mark.parametrize("material", ["matte", "plastic", "glossy", "metallic"])
def test_material_is_camera_relative_and_two_sided(material):
    p = profile(material)
    color = [[0.1, 0.6, 0.45]]
    front = shade(color, [[0, 0, 1]], [[0, 0, 1]], p)
    back = shade(color, [[0, 0, -1]], [[0, 0, 1]], p)
    np.testing.assert_allclose(front, back)
    assert np.isfinite(front).all()


def test_sheet_smoothing_uses_neighbors_before_interpolation():
    points = np.array([[0, 0, 0], [3, 1, 0], [6, -1, 0], [9, 0, 0]], float)
    control = controls_for(points, np.tile([0, 0, 1], (4, 1)), ["", "S", "S", "H"], 1)
    np.testing.assert_allclose(control["p2"], [3, 0.25, 0])
    assert control["secStrucFirst"]
    assert not control["secStrucLast"]


def test_residue_curves_share_boundary_positions_and_frames():
    points = np.array([[0, 0, 0], [3, 1, 1], [6, -1, 1], [9, 0, 0], [12, 1, 0]], float)
    directions = np.tile([0, 0, 1], (5, 1))
    a = interpolate(controls_for(points, directions, ["H"] * 5, 1), 0.9)
    b = interpolate(controls_for(points, directions, ["H"] * 5, 2), 0.9)
    for first, second in zip(a, b):
        np.testing.assert_allclose(first[-1], second[0], atol=1e-6)


def test_beta_cross_section_has_flat_faces_and_sharp_corners():
    curve = np.array([[0, 0, 0], [2, 0, 0]], float)
    normals = np.tile([0, 1, 0], (2, 1))
    binormals = np.tile([0, 0, 1], (2, 1))
    mesh = section_mesh(
        curve,
        normals,
        binormals,
        np.full(2, 0.2),
        np.ones(2),
        [0.1, 0.6, 0.4],
        7,
        "square",
        16,
        (False, False),
    )
    np.testing.assert_allclose(abs(mesh.vertices[:, 1]), 1)
    np.testing.assert_allclose(abs(mesh.vertices[:, 2]), 0.2)
    assert len(np.unique(mesh.normals, axis=0)) == 4
    assert np.all(mesh.owners == 7)


def test_elliptical_helix_and_square_sheet_are_distinct_profiles():
    curve = np.array([[0, 0, 0], [2, 0, 0]], float)
    n = np.tile([0, 1, 0], (2, 1))
    b = np.tile([0, 0, 1], (2, 1))
    mesh = section_mesh(
        curve,
        n,
        b,
        np.full(2, 0.2),
        np.ones(2),
        [1, 0, 0],
        0,
        "elliptical",
        16,
        (False, False),
    )
    np.testing.assert_allclose(
        mesh.vertices[:, 1] ** 2 + (mesh.vertices[:, 2] / 0.2) ** 2, 1, atol=1e-6
    )


def test_ssao_samples_are_deterministic_hemisphere_vectors():
    samples = ssao_samples()
    assert samples.shape == (32, 3)
    assert np.all(samples[:, 2] >= 0)
    assert np.all(
        np.linalg.norm(samples, axis=1)
        <= 0.1 + 0.9 * ((np.arange(32) + 1) / 32) ** 2 + 1e-6
    )
    np.testing.assert_array_equal(samples, ssao_samples())


def test_flat_plane_has_no_spurious_occlusion_in_headless_export():
    from molstar_style_in_pymol.occlusion import calculate

    cmd = SimpleNamespace(
        get_view=lambda: [*np.eye(3).ravel(), 0, 0, -20, 0, 0, 0, 1, 40, 45],
        get_setting_int=lambda _: 1,
        get_setting_float=lambda _: 45,
    )
    distance = np.full((32, 48), 20.0)
    np.testing.assert_allclose(calculate(distance, cmd, {}), 1, atol=1e-6)
    distance[10:22, 15:30] = 17
    result = calculate(distance, cmd, {})
    assert np.isfinite(result).all()
    assert result.min() < 0.98
    assert np.all((result >= 0.01) & (result <= 1))


def test_mixed_secondary_structure_keeps_sheet_faces_flat():
    from dataclasses import replace

    from test_reference_geometry import atoms_at, draw

    state = atoms_at([[3.8 * i, 0, 0] for i in range(5)], kind="protein")
    state.atoms = tuple(
        replace(a, ss=ss) for a, ss in zip(state.atoms, ["H", "S", "S", "S", "H"])
    )
    mesh = (
        draw(state, "cartoon", {"radialSegments": 16, "linearSegments": 4})
        .pieces[0]
        .mesh
    )
    normals = mesh.normals[mesh.owners == 2]
    assert len(np.unique(np.round(normals, 5), axis=0)) == 4
