"""Scientific value/format checks independent of renderer registration."""

import struct

import numpy as np
import pytest
from samples import molecule, scalar_grid

from molstar_style_in_pymol.data import atom_fields, local_path, read_cif, read_grid
from molstar_style_in_pymol.registry import reference
from molstar_style_in_pymol.themes import ANNOTATION_THEMES, colors, sizes


def annotation_data(state):
    n = len(state.atoms)
    keys = (
        "entity_id",
        "entity_source",
        "operator_hkl",
        "operator_name",
        "group",
        "instance",
        "segment",
        "compartment",
        "entity",
        "hierarchy",
        "particle_index",
        "value",
        "radius",
        "size",
        *(k.replace("-", "_") for k in ANNOTATION_THEMES),
    )
    return {
        "atom_data": {k: np.arange(n) % 5 + 1 for k in keys},
        "grid": scalar_grid(),
        "external_coordinates": state.coords,
        "external_colors": [a.color for a in state.atoms],
        "layers": [{"field": "value", "color": "red"}],
    }


@pytest.mark.parametrize(
    "theme", (*reference()["color_themes"], *ANNOTATION_THEMES, "auto", "keep")
)
def test_color_themes(theme):
    state = molecule()
    c = colors(state.atoms, state.coords, theme, data=annotation_data(state))
    assert c.shape == (len(state.atoms), 3)
    assert np.isfinite(c).all() and c.min() >= 0 and c.max() <= 1


@pytest.mark.parametrize("theme", reference()["size_themes"])
def test_size_themes(theme):
    state = molecule()
    r = sizes(state.atoms, theme, data=annotation_data(state))
    assert len(r) == len(state.atoms) and np.isfinite(r).all() and r.min() >= 0


def test_no_runtime_network():
    for url in (
        "https://example.org/1abc.pdb",
        "http://example.org/grid",
        "ftp://example.org/file",
    ):
        with pytest.raises(ValueError, match="local"):
            local_path(url)


def test_map_ccp4_axis_and_origin(tmp_path):
    header = bytearray(1024)
    for index, value in {
        0: 3,
        1: 4,
        2: 5,
        3: 2,
        4: 1,
        5: 2,
        6: 3,
        7: 4,
        8: 3,
        9: 5,
        16: 2,
        17: 1,
        18: 3,
        22: 1,
    }.items():
        struct.pack_into("<i", header, index * 4, value)
    for index, value in {10: 8, 11: 9, 12: 10, 13: 90, 14: 90, 15: 90}.items():
        struct.pack_into("<f", header, index * 4, value)
    header[208:212] = b"MAP "
    values = np.arange(60, dtype=np.float32).reshape(5, 4, 3)
    path = tmp_path / "density.mrc"
    path.write_bytes(header + values.tobytes())
    grid = read_grid(path)
    np.testing.assert_array_equal(grid.values, values.transpose(2, 1, 0))
    np.testing.assert_allclose(
        grid.transform[:3, :3], [[0, 2, 0], [3, 0, 0], [0, 0, 2]], atol=1e-6
    )
    np.testing.assert_allclose(grid.world([[0, 0, 0]])[0], [4, 3, 6])


def test_cube_negative_axis_units(tmp_path):
    path = tmp_path / "density.cube"
    path.write_text(
        "comment\ncomment\n0 1 2 3\n-2 0.5 0 0\n-2 0 0.6 0\n-2 0 0 0.7\n0 1 2 3 4 5 6 7\n"
    )
    grid = read_grid(path)
    np.testing.assert_allclose(grid.spacing, [0.5, 0.6, 0.7])
    assert grid.values[1, 1, 1] == 7


def test_bcif_fixed_point_and_mask(tmp_path):
    import msgpack

    payload = {
        "dataBlocks": [
            {
                "categories": [
                    {
                        "name": "_test",
                        "rowCount": 3,
                        "columns": [
                            {
                                "name": "value",
                                "data": {
                                    "data": np.array([10, 20, 30], "<i4").tobytes(),
                                    "encoding": [
                                        {"kind": "FixedPoint", "factor": 10},
                                        {"kind": "ByteArray", "type": 3},
                                    ],
                                },
                                "mask": {
                                    "data": bytes([0, 1, 0]),
                                    "encoding": [{"kind": "ByteArray", "type": 4}],
                                },
                            }
                        ],
                    }
                ]
            }
        ]
    }
    path = tmp_path / "sample.bcif"
    path.write_bytes(msgpack.packb(payload))
    rows = read_cif(path)["test"]
    assert rows == [{"value": 1.0}, {"value": None}, {"value": 3.0}]


def test_duplicate_residue_annotation_rejected():
    state = molecule()
    row = {"chain": "A", "resi": "1", "plddt": 90}
    with pytest.raises(ValueError, match="Ambiguous"):
        atom_fields(state.atoms, {"residues": [row, row]})


def test_mvs_local_structure(tmp_path):
    import pymol2

    from molstar_style_in_pymol.geometry import build

    path = tmp_path / "model.pdb"
    state = molecule()
    with pymol2.PyMOL() as instance:
        instance.cmd.load_model(state.model, "sample")
        instance.cmd.save(str(path), "sample")
    data = {
        "_base": str(tmp_path),
        "root": {
            "kind": "root",
            "children": [
                {
                    "kind": "download",
                    "params": {"url": "model.pdb"},
                    "children": [
                        {
                            "kind": "parse",
                            "params": {"format": "pdb"},
                            "children": [
                                {
                                    "kind": "structure",
                                    "params": {"type": "model"},
                                    "children": [
                                        {
                                            "kind": "component",
                                            "params": {"selector": "all"},
                                            "children": [
                                                {
                                                    "kind": "representation",
                                                    "params": {
                                                        "type": "ball-and-stick"
                                                    },
                                                    "children": [
                                                        {
                                                            "kind": "color",
                                                            "params": {"color": "red"},
                                                        }
                                                    ],
                                                }
                                            ],
                                        }
                                    ],
                                }
                            ],
                        }
                    ],
                }
            ],
        },
    }
    g = build(None, "mvs", "auto", {}, data, "low")
    assert g.pieces and g.nbytes > 1000
    assert g.pieces[0].mesh.colors[:, 0].min() == 1


def test_saccharide_reference_size():
    assert len(reference()["saccharides"]) > 2000
    assert reference()["saccharides"]["NAG"]["shape"] == "FilledCube"
