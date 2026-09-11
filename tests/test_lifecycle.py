"""Real PyMOL transaction, overlapping views, and session round trips."""

from unittest.mock import patch

import numpy as np
import pymol2
import pytest
from samples import molecule

from molstar_style_in_pymol import molstar_style
from molstar_style_in_pymol.controller import manager_for


@pytest.fixture
def cmd():
    with pymol2.PyMOL() as instance:
        instance.cmd.load_model(molecule().model, "sample")
        instance.cmd.show_as("lines", "sample")
        yield instance.cmd


def reps(cmd):
    rows = []
    cmd.iterate("sample", "rows.append((index,reps,color,ss))", space={"rows": rows})
    return rows


def test_overlap_and_reset(cmd):
    before = reps(cmd)
    coords = cmd.get_coords("sample").copy()
    molstar_style("spacefill", "sample", name="view_a", _self=cmd)
    molstar_style("ball-and-stick", "sample and chain A", name="view_b", _self=cmd)
    molstar_style("reset", name="view_a", _self=cmd)
    assert all(row[1] == 0 for row in reps(cmd)[:24])
    molstar_style("reset", name="view_b", _self=cmd)
    assert reps(cmd) == before
    np.testing.assert_array_equal(cmd.get_coords("sample"), coords)


def test_failed_preparation_keeps_existing_view(cmd):
    entry = molstar_style("spacefill", _self=cmd)
    objects = cmd.get_names("objects")
    with pytest.raises(Exception, match="requires"):
        molstar_style("quality-plddt", _self=cmd)
    assert manager_for(cmd).entries["molstar"] is entry
    assert cmd.get_names("objects") == objects


def test_load_failure_rolls_back(cmd):
    before = reps(cmd)
    with patch.object(
        cmd, "load_cgo", side_effect=RuntimeError("injected load failure")
    ):
        with pytest.raises(Exception, match="injected load failure"):
            molstar_style("spacefill", _self=cmd)
    assert reps(cmd) == before
    assert cmd.get_names("objects") == ["sample"]


def test_manual_delete_restores(cmd):
    before = reps(cmd)
    entry = molstar_style("cartoon", _self=cmd)
    cmd.delete(entry.generated[0])
    manager_for(cmd).maintenance()
    assert not manager_for(cmd).entries
    assert reps(cmd) == before


def test_session_round_trip(cmd):
    before = reps(cmd)
    molstar_style("spacefill", _self=cmd)
    session = cmd.get_session()
    cmd.reinitialize()
    cmd.set_session(session)
    assert "molstar" in manager_for(cmd).entries
    molstar_style("reset", _self=cmd)
    assert reps(cmd) == before


def test_all_states_and_refresh(cmd):
    xyz = cmd.get_coords("sample")
    cmd.load_coordset(xyz + [0, 0, 2], "sample", state=2)
    entry = molstar_style("cartoon", _self=cmd)
    assert len(entry.drawings["selection"]) == 2
    expected = np.array(cmd.get_model("sample", state=2).get_coord_list())
    np.testing.assert_allclose(entry.drawings["selection"][1].anchors, expected)
    cmd.translate([1, 0, 0], "sample", state=1, camera=0)
    molstar_style("refresh", _self=cmd)
    expected = np.array(cmd.get_model("sample", state=1).get_coord_list())
    np.testing.assert_allclose(
        manager_for(cmd).entries["molstar"].drawings["selection"][0].anchors, expected
    )


def test_import_has_no_pymol_side_effect():
    import subprocess
    import sys

    subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import molstar_style_in_pymol; assert 'pymol' not in sys.modules",
        ],
        check=True,
    )


def test_volume_ray_accumulates_opacity_and_restores_settings(cmd, tmp_path):
    from PIL import Image
    from samples import scalar_grid

    settings = (
        "transparency_mode",
        "triangle_max_passes",
        "ray_max_passes",
        "backface_cull",
    )
    before = {k: cmd.get_setting_tuple(k) for k in settings}
    cmd.hide("everything", "sample")
    entry = molstar_style(
        "direct-volume", "none", data={"grid": scalar_grid()}, _self=cmd
    )
    cmd.orient(entry.name)
    cmd.zoom(entry.name, 1)
    path = tmp_path / "density.png"
    molstar_style("ray", filename=str(path), width=120, height=100, _self=cmd)
    image = Image.open(path).convert("RGBA")
    assert np.asarray(image)[:, :, 3].max() > 200
    assert {k: cmd.get_setting_tuple(k) for k in settings} == before
    cmd.ray(120, 100, quiet=1)
    from io import BytesIO

    native = Image.open(BytesIO(cmd.png(None, prior=1, quiet=1))).convert("RGBA")
    assert np.asarray(native)[:, :, 3].max() > 80


def test_source_state_and_visibility_follow(cmd):
    cmd.load_coordset(cmd.get_coords("sample") + [0, 0, 1], "sample", state=2)
    entry = molstar_style("cartoon", _self=cmd)
    cmd.set("state", 2, "sample")
    manager_for(cmd).maintenance()
    assert all(cmd.get_setting_int("state", n) == 2 for n in entry.generated)
    cmd.disable("sample")
    manager_for(cmd).maintenance()
    assert not list(manager_for(cmd).active_drawings())
    cmd.enable("sample")
    manager_for(cmd).maintenance()
    assert list(manager_for(cmd).active_drawings())
