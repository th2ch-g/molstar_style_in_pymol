"""Analytic integration and occlusion checks for headless density rendering."""

from types import SimpleNamespace

import numpy as np

from molstar_style_in_pymol.data import Grid
from molstar_style_in_pymol.scene import Volume
from molstar_style_in_pymol.volume_ray import render


def scene():
    view = [*np.eye(3).ravel(), 0, 0, -10, 0, 0, 0, 1, 20, 20]
    cmd = SimpleNamespace(
        get_view=lambda: view,
        get_setting_int=lambda key: 1,
        get_setting_float=lambda key: 20,
    )
    transform = np.eye(4)
    transform[:3, 3] = [-1, -1, -2]
    volume = Volume(
        Grid(np.full((3, 3, 5), 0.5), transform),
        np.array([[0, 0.1, 0.6, 0.4, 0.2], [1, 0.1, 0.6, 0.4, 0.2]]),
        step=0.5,
    )
    return cmd, volume


def test_constant_density_matches_beer_lambert_sampling():
    cmd, volume = scene()
    pixel = render(volume, cmd, np.full((3, 3), np.inf), 3, 3)[1, 1]
    alpha = 1 - 0.8**9
    np.testing.assert_allclose(
        pixel, [0.1 * alpha, 0.6 * alpha, 0.4 * alpha, alpha], atol=1e-7
    )


def test_native_opaque_depth_clips_density_integration():
    cmd, volume = scene()
    hidden = render(volume, cmd, np.full((3, 3), 7), 3, 3)
    np.testing.assert_array_equal(hidden, 0)
    partial = render(volume, cmd, np.full((3, 3), 9), 3, 3)[1, 1]
    np.testing.assert_allclose(partial[3], 1 - 0.8**3, atol=1e-7)
