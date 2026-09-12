"""Numerical regressions against the pinned Mol* formulas and geometric invariants."""

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest
from samples import data_for, scalar_grid
from scipy.spatial.transform import Rotation

from molstar_style_in_pymol.data import Grid
from molstar_style_in_pymol.geometry import build
from molstar_style_in_pymol.molecule import axes_of, ellipsoid
from molstar_style_in_pymol.source import Atom, State
from molstar_style_in_pymol.volume import density_projection, density_slices, isosurface


def atoms_at(
    coords, bonds=(), aniso=(0.04, 0.09, 0.16, 0.02, 0.01, 0.015), kind="other", order=1
):
    atoms = tuple(
        Atom(
            "model",
            i + 1,
            "CA",
            str(i + 1),
            "ALA",
            "A",
            "",
            "C",
            "",
            kind,
            (1, 0, 0) if i % 2 == 0 else (0, 0, 1),
            1.7,
            bfactor=25,
            aniso=aniso,
        )
        for i in range(len(coords))
    )
    model = SimpleNamespace(
        bond=[SimpleNamespace(order=order, index=list(b)) for b in bonds]
    )
    return State(
        atoms,
        np.array(coords, float),
        tuple(bonds),
        model,
        np.zeros(len(coords), int),
        {},
    )


def draw(state, rep, params=None, data=None, color="auto"):
    return build(state, rep, color, params or {}, data or {}, "high")


def face_alignment(mesh):
    triangles = mesh.vertices[mesh.faces]
    face_normals = np.cross(
        triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
    )
    return np.einsum("ij,ij->i", face_normals, mesh.normals[mesh.faces].mean(axis=1))


def test_anisotropic_ellipsoid_equation_and_normals():
    state = atoms_at([[3, -2, 4]])
    mesh = draw(state, "ellipsoid").pieces[0].mesh
    u11, u22, u33, u12, u13, u23 = state.atoms[0].aniso
    tensor = np.array([[u11, u12, u13], [u12, u22, u23], [u13, u23, u33]])
    delta = mesh.vertices - state.coords[0]
    np.testing.assert_allclose(
        np.einsum("ij,jk,ik->i", delta, np.linalg.inv(tensor), delta),
        1.5958**2,
        rtol=2e-6,
    )
    assert face_alignment(mesh).min() > 0
    analytical = delta @ np.linalg.inv(tensor)
    analytical /= np.linalg.norm(analytical, axis=1)[:, None]
    np.testing.assert_allclose(mesh.normals, analytical, atol=2e-6)


@pytest.mark.parametrize(
    "params,radius",
    [
        ({}, 0.1),
        ({"sizeFactor": 2}, 0.2),
        ({"sizeTheme": "physical"}, 0.17),
        ({"sizeAspectRatio": 0.2}, 0.2),
    ],
)
def test_ellipsoid_bond_radius(params, radius):
    state = atoms_at([[0, 0, 0], [3, 0, 0]], [(0, 1)])
    mesh = (
        draw(state, "ellipsoid", {**params, "visuals": ["intra-bond"]}).pieces[0].mesh
    )
    np.testing.assert_allclose(
        np.linalg.norm(mesh.vertices[:, 1:], axis=1).max(), radius, rtol=1e-6
    )


def test_ellipsoid_missing_and_isotropic_tensors():
    state = atoms_at([[0, 0, 0], [5, 0, 0]], aniso=(0.09, 0.09, 0.09, 0, 0, 0))
    state.atoms = (state.atoms[0], replace(state.atoms[1], aniso=()))
    mesh = draw(state, "ellipsoid").pieces[0].mesh
    assert set(mesh.owners) == {0}
    np.testing.assert_allclose(
        np.linalg.norm(mesh.vertices, axis=1), 0.3 * 1.5958, atol=1e-6
    )
    with pytest.raises(ValueError, match="anisotropic"):
        draw(atoms_at([[0, 0, 0]], aniso=()), "ellipsoid")


def test_explicit_ellipsoid_probability():
    from scipy.stats import chi2

    mesh = (
        draw(
            atoms_at([[0, 0, 0]], aniso=(1, 1, 1, 0, 0, 0)),
            "ellipsoid",
            {"probability": 0.5},
        )
        .pieces[0]
        .mesh
    )
    np.testing.assert_allclose(
        np.linalg.norm(mesh.vertices, axis=1), np.sqrt(chi2.ppf(0.5, 3)), atol=1e-6
    )


def test_reflected_ellipsoid_keeps_outward_winding():
    normal = ellipsoid([0, 0, 0], np.eye(3), [1, 2, 3], [1, 0, 0], 0, 12)
    original = normal.faces.copy()
    reflected = ellipsoid([0, 0, 0], np.diag([-1, 1, 1]), [1, 2, 3], [1, 0, 0], 0, 12)
    assert face_alignment(reflected).min() > 0
    np.testing.assert_array_equal(normal.faces, original)


@pytest.mark.parametrize(
    "rep,radius", [("backbone", 0.3), ("putty", 0.54), ("cartoon", 0.2)]
)
def test_polymer_default_radius(rep, radius):
    mesh = (
        draw(atoms_at([[0, 0, 0], [3.8, 0, 0], [7.6, 0, 0]], kind="protein"), rep)
        .pieces[0]
        .mesh
    )
    np.testing.assert_allclose(
        np.linalg.norm(mesh.vertices[:, 1:], axis=1).max(), radius, atol=1e-6
    )


@pytest.mark.parametrize("rep", ["backbone", "cartoon", "putty"])
def test_polymer_size_theme_is_applied(rep):
    state = atoms_at([[0, 0, 0], [3.8, 0, 0]], kind="protein")
    params = {"sizeTheme": "uniform", "sizeParams": {"value": 3}}
    mesh = draw(state, rep, params).pieces[0].mesh
    expected = 0.9 if rep == "backbone" else 0.6
    np.testing.assert_allclose(
        np.linalg.norm(mesh.vertices[:, 1:], axis=1).max(), expected, atol=1e-6
    )


def test_line_attenuation_boolean_is_not_a_radius():
    state = atoms_at([[0, 0, 0], [3, 0, 0]], [(0, 1)])
    for flag in (False, True):
        mesh = (
            draw(
                state, "line", {"visuals": ["intra-bond"], "lineSizeAttenuation": flag}
            )
            .pieces[0]
            .mesh
        )
        np.testing.assert_allclose(
            np.linalg.norm(mesh.vertices[:, 1:], axis=1).max(), 0.04, atol=1e-6
        )
    assert draw(atoms_at([[0, 0, 0]]), "line").pieces
    assert not draw(state, "line", {"visuals": ["element-cross"]}).pieces
    assert draw(state, "line", {"visuals": ["element-cross"], "crosses": "all"}).pieces


def test_double_bond_radius_and_plane():
    state = atoms_at([[0, 0, 0], [3, 0, 0], [0, 2, 0]], [(0, 1), (0, 2)], order=2)
    mesh = draw(state, "ball-and-stick", {"visuals": ["intra-bond"]}).pieces[0].mesh
    np.testing.assert_allclose(
        abs(mesh.vertices[:, 2]).max(), 1.7 * 0.15 * (2 / 3) * 0.45, atol=1e-6
    )


@pytest.mark.parametrize("sign", [1, -1])
@pytest.mark.parametrize("reflection", [1, -1])
def test_surface_winding_and_negative_lobes(sign, reflection):
    indices = np.indices((17, 17, 17))
    field = sign * np.exp(-((indices - 8) ** 2).sum(axis=0) / 12)
    transform = np.diag([reflection * 0.4, 0.7, 1.1, 1])
    grid = Grid(field, transform)
    mesh = isosurface(grid, sign * 0.5)
    assert face_alignment(mesh).min() > 0
    outward = mesh.vertices - grid.world([[8, 8, 8]])[0]
    assert np.einsum("ij,ij->i", outward, mesh.normals).min() > 0


def test_pca_box_midpoint_and_transformation():
    points = np.array([[0, 0, 0], [0.2, 0.1, 0], [0.3, -0.1, 0.1], [6, 2, 1]])
    for rotation in (
        np.eye(3),
        Rotation.from_euler("xyz", [20, -35, 12], degrees=True).as_matrix(),
    ):
        transformed = points @ rotation.T + [5, -8, 3]
        center, axes, extent = axes_of(transformed)
        projected = (transformed - center) @ axes.T
        np.testing.assert_allclose(
            projected.max(axis=0) + projected.min(axis=0), 0, atol=1e-10
        )
        assert (abs(projected) <= extent + 1e-10).all()
        assert np.linalg.det(axes) > 0


def test_orientation_default_only_ellipsoid_and_scaling():
    state = atoms_at([[0, 0, 0], [3, 2, 1], [1, 4, -1]])
    geometry = draw(state, "orientation")
    assert len(geometry.pieces) == 1
    center, _, _ = axes_of(state.coords)
    large = draw(state, "orientation", {"sizeFactor": 2}).pieces[0].mesh
    np.testing.assert_allclose(
        large.vertices - center,
        2 * (geometry.pieces[0].mesh.vertices - center),
        atol=1e-6,
    )
    default = (
        draw(None, "shape-orientation", data={"positions": state.coords}).pieces[0].mesh
    )
    box = (
        draw(
            None, "shape-orientation", {"visuals": ["box"]}, {"positions": state.coords}
        )
        .pieces[0]
        .mesh
    )
    np.testing.assert_array_equal(default.vertices, box.vertices)


def test_structure_plane_cross_section_and_cutout():
    state = atoms_at([[0, 0, 0]])
    params = {
        "mode": "plane",
        "plane": {"point": [0, 0, 1], "normal": [0, 0, 1]},
        "cutout": True,
        "antialias": False,
        "imageResolution": 0.04,
    }
    mesh = draw(state, "plane", params, color="keep").pieces[0].mesh
    np.testing.assert_allclose(mesh.vertices[:, 2], 1, atol=1e-6)
    radius = np.linalg.norm(mesh.vertices[:, :2], axis=1).max()
    assert abs(radius - np.sqrt(1.7**2 - 1)) < 0.05
    np.testing.assert_allclose(mesh.colors, np.tile([1, 0, 0], (len(mesh.vertices), 1)))
    params["plane"]["point"] = [0, 0, 3]
    assert not draw(state, "plane", params).pieces
    params["cutout"] = False
    assert draw(state, "plane", params).pieces


@pytest.mark.parametrize(
    "rep", ["dot", "slice", "isosurface", "segment", "direct-volume"]
)
def test_volume_uniform_color(rep):
    geometry = draw(
        None,
        rep,
        {"colorParams": {"value": 0xFF0000}, "sizeFactor": 0.1},
        data_for(rep),
        "uniform",
    )
    for piece in geometry.pieces:
        if rep == "slice":
            np.testing.assert_array_equal(piece.mesh.colors[:, 1:], 0)
            assert np.ptp(piece.mesh.colors[:, 0]) > 0.5
        else:
            np.testing.assert_allclose(
                piece.mesh.colors, np.tile([1, 0, 0], (len(piece.mesh.vertices), 1))
            )
    for volume in geometry.volumes:
        np.testing.assert_allclose(
            volume.lookup[:, :3], np.tile([1, 0, 0], (len(volume.lookup), 1))
        )


def test_gaussian_volume_atom_colors_survive_both_ray_paths():
    state = atoms_at([[0, 0, 0], [6, 0, 0]])
    volume = draw(state, "gaussian-volume", {"resolution": 0.5}, color="keep").volumes[
        0
    ]
    np.testing.assert_allclose(
        volume.sample_colors(state.coords), [[1, 0, 0], [0, 0, 1]]
    )
    assert volume.pixels.shape == (*volume.grid.values.shape[::-1], 4)
    for mesh in (
        density_projection(volume),
        density_slices(volume, direction=[0, 0, 1]),
    ):
        assert mesh.colors[:, 0].max() > 0.9 and mesh.colors[:, 2].max() > 0.9
        assert mesh.colors[:, 1].max() == 0


def test_negative_volume_dot_threshold():
    grid = Grid(
        np.array([[[-2, -1, 0], [-1, 0, 1]], [[0, 0, 0], [0, 0, 0]]]), np.eye(4)
    )
    mesh = (
        draw(None, "dot", {"isoValue": -1.5, "sizeFactor": 0.1}, {"grid": grid})
        .pieces[0]
        .mesh
    )
    assert np.linalg.norm(mesh.vertices, axis=1).max() < 0.101


@pytest.mark.parametrize("degrees", [0, 90, 180])
def test_angle_sector_and_collinear_limits(degrees):
    radians = np.deg2rad(degrees)
    positions = [[2, 0, 0], [0, 0, 0], [2 * np.cos(radians), 2 * np.sin(radians), 0]]
    geometry = draw(
        None, "angle", {"visuals": ["sector"], "label": False}, {"positions": positions}
    )
    assert geometry.metadata["value"] == pytest.approx(degrees)
    if degrees:
        mesh = geometry.pieces[0].mesh
        np.testing.assert_allclose(
            np.linalg.norm(mesh.vertices, axis=1).max(), 1.4, atol=1e-6
        )
        triangles = mesh.vertices[mesh.faces]
        area = (
            np.linalg.norm(
                np.cross(
                    triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]
                ),
                axis=1,
            ).sum()
            / 2
        )
        assert area == pytest.approx(0.5 * radians * 1.4**2, rel=0.002)


def test_confal_midpoint_and_reference_colors():
    positions = np.array([[0, 1, 2], [0, 0, 0], [-1, 0, 0], [1, 0, 0], [0, -1, 2]])
    data = {"steps": [{"positions": positions, "class": "AB01", "score": 20}]}
    mesh = draw(None, "confal-pyramids", data=data).pieces[0].mesh
    assert np.any(np.all(np.isclose(mesh.vertices, [0, 0, 2]), axis=1))
    assert not np.any(np.all(np.isclose(mesh.vertices, [0, 0, 0]), axis=1))
    colors = np.unique(mesh.colors, axis=0)
    assert len(colors) == 2
    assert any(np.allclose(c, np.array([255, 193, 193]) / 255) for c in colors)


def test_mesh_cif_noncontiguous_vertex_ids():
    data = {
        "tables": {
            "mesh": [{"id": 7}],
            "mesh_vertex": [
                {"mesh_id": 7, "vertex_id": key, "x": pos[0], "y": pos[1], "z": pos[2]}
                for key, pos in ((42, [0, 0, 0]), (9, [1, 0, 0]), (100, [0, 1, 0]))
            ],
            "mesh_triangle": [{"mesh_id": 7, "vertex_id": key} for key in (42, 9, 100)],
        }
    }
    mesh = draw(None, "mesh", data=data).pieces[0].mesh
    np.testing.assert_array_equal(mesh.faces, [[0, 1, 2]])


def test_particle_target_instances_rotate_center_and_scale():
    data = data_for("particle-target")
    target = data["targets"]["sample"]
    target["vertices"] = (np.array(target["vertices"]) + [10, 20, 30]).tolist()
    quaternion = Rotation.from_euler("z", 90, degrees=True).as_quat()
    data["particles"] = [
        {
            "target": "sample",
            "position": [3, 5, 7],
            "radius": 2,
            "quaternion": quaternion,
        }
    ]
    mesh = draw(None, "particle-target", data=data).pieces[0].mesh
    expected = (np.array(target["vertices"]) - [10, 20, 30]) @ Rotation.from_quat(
        quaternion
    ).as_matrix().T * 2 + [3, 5, 7]
    np.testing.assert_allclose(mesh.vertices, expected, atol=1e-6)
    assert len(mesh.faces) == 4


@pytest.mark.parametrize("kind", ["structure", "volume"])
def test_particle_structure_and_volume_targets(kind):
    target = {
        "kind": kind,
        "positions": [[0, 0, 0]],
        "radii": [1],
        "grid": scalar_grid(),
    }
    data = {
        "targets": {"a": target},
        "particles": [{"target": "a", "position": [10, 0, 0]}],
    }
    mesh = draw(None, "particle-target", data=data).pieces[0].mesh
    np.testing.assert_allclose(
        (mesh.vertices.min(axis=0) + mesh.vertices.max(axis=0)) / 2,
        [10, 0, 0],
        atol=1e-5,
    )


def test_particle_scale_and_orientation_length():
    data = {"particles": [{"position": [0, 0, 0], "radius": 2, "scale": [1, 2, 3]}]}
    mesh = draw(None, "particle-spacefill", {"sizeFactor": 2}, data).pieces[0].mesh
    np.testing.assert_allclose(abs(mesh.vertices).max(axis=0), [4, 8, 12], atol=1e-6)
    mesh = draw(None, "particle-orientation", data=data).pieces[0].mesh
    np.testing.assert_allclose(mesh.vertices.max(axis=0), [10, 10, 10], atol=1e-6)


def test_triclinic_unitcell_volume():
    data = {"cell": [10, 12, 14, 90, 100, 110]}
    mesh = draw(None, "unitcell", data=data).pieces[0].mesh
    assert np.isfinite(mesh.vertices).all()
    for cell in ([0, 12, 14, 90, 90, 90], [10, 12, 14, 0, 90, 90]):
        with pytest.raises(ValueError):
            draw(None, "unitcell", data={"cell": cell})


def test_ray_depth_rounding_at_image_border():
    from unittest.mock import patch

    from molstar_style_in_pymol.ray_effects import depth_samples

    view = [0] * 18
    view[11] = -10
    cmd = SimpleNamespace(
        get_view=lambda: view,
        get_setting_float=lambda _: 90,
        get_setting_int=lambda _: 1,
    )
    mesh = SimpleNamespace(
        vertices=np.array([[-1, -10.07, -10], [1, -10.07, -10], [0, -9, -9]]),
        faces=np.array([[0, 1, 2]]),
    )
    drawing = SimpleNamespace(pieces=[SimpleNamespace(mesh=mesh)])
    with patch("molstar_style_in_pymol.export.view_matrix", return_value=np.eye(4)):
        depth = depth_samples([drawing], cmd, 1200, 900)
    assert depth.shape == (900, 1200) and np.isfinite(depth).all()


def test_ray_depth_raster_fills_triangle_interiors():
    from unittest.mock import patch

    from molstar_style_in_pymol.ray_effects import depth_samples

    view = [0] * 18
    view[11] = -10
    cmd = SimpleNamespace(
        get_view=lambda: view,
        get_setting_float=lambda _: 90,
        get_setting_int=lambda _: 1,
    )
    mesh = SimpleNamespace(
        vertices=np.array([[-5, -5, -10], [5, -5, -10], [5, 5, -10], [-5, 5, -10]]),
        faces=np.array([[0, 1, 2], [0, 2, 3]]),
    )
    with patch("molstar_style_in_pymol.export.view_matrix", return_value=np.eye(4)):
        depth = depth_samples(
            [SimpleNamespace(pieces=[SimpleNamespace(mesh=mesh)])], cmd, 80, 60
        )
    np.testing.assert_allclose(depth[17:43, 28:52], 0.1, atol=1e-6)
    np.testing.assert_allclose(depth[:10], 1)


def test_ray_occlusion_preserves_flat_background():
    from io import BytesIO
    from unittest.mock import patch

    from PIL import Image

    from molstar_style_in_pymol.ray_effects import process

    pixels = np.full((60, 80, 4), 255, np.uint8)
    pixels[20:40, 20:60, 1:3] = 0
    stream = BytesIO()
    Image.fromarray(pixels).save(stream, format="PNG")
    cmd = SimpleNamespace(
        get_setting_int=lambda _: 0, get_setting_tuple=lambda _: (4, (1, 1, 1))
    )
    profile = SimpleNamespace(effects={"occlusion": True}, background=None)
    depth = np.full((60, 80), 0.2)
    depth[20:40, 20:60] = 0.1
    with patch("molstar_style_in_pymol.ray_effects.depth_samples", return_value=depth):
        output = process(
            stream.getvalue(),
            [SimpleNamespace(profile=profile, volumes=[])],
            cmd,
            ambient_occlusion=depth,
        )
    image = np.asarray(Image.open(BytesIO(output)))
    background = np.all(pixels[:, :, :3] == 255, axis=2)
    np.testing.assert_array_equal(image[background], pixels[background])
