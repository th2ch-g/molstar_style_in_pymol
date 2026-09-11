"""Local scientific inputs and portable, non-pickle scene recipes."""

import gzip
import itertools
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


def local_path(value, base=None):
    text = str(value)
    if re.match(r"^[a-zA-Z][a-zA-Z0-9+.-]*://", text):
        raise ValueError(
            "Only local files are supported; download remote assets separately"
        )
    path = Path(text).expanduser()
    if base is not None and not path.is_absolute():
        path = Path(base) / path
    if not path.is_file():
        raise ValueError(f"Input file does not exist: {path.name}")
    return path


def parameters(value):
    if value is None or value == "":
        return {}
    if isinstance(value, dict):
        return json.loads(json.dumps(value, allow_nan=False))
    path = local_path(value)
    result = json.loads(path.read_text())
    if not isinstance(result, dict):
        raise ValueError("params must be a JSON object")
    return result


def finite(value, shape=None, name="array", dtype=float):
    array = np.asarray(value, dtype=dtype)
    if shape is not None and (
        array.ndim != len(shape)
        or any(n is not None and n != a for n, a in zip(shape, array.shape))
    ):
        raise ValueError(f"{name} must have shape {shape}; got {array.shape}")
    if not np.isfinite(array).all():
        raise ValueError(f"{name} contains non-finite values")
    return array


@dataclass
class Grid:
    values: np.ndarray
    transform: np.ndarray
    label: str = "grid"
    metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        self.values = finite(self.values, (None, None, None), "grid", np.float32)
        if min(self.values.shape) < 1:
            raise ValueError("Grid dimensions must be nonzero")
        self.transform = finite(self.transform, (4, 4), "grid transform")
        if abs(np.linalg.det(self.transform[:3, :3])) < 1e-12:
            raise ValueError("Grid transform is singular")

    @property
    def nbytes(self):
        return self.values.nbytes + self.transform.nbytes

    @property
    def spacing(self):
        return np.linalg.norm(self.transform[:3, :3], axis=0)

    def world(self, points):
        return np.asarray(points) @ self.transform[:3, :3].T + self.transform[:3, 3]

    def sample(self, points):
        from scipy.ndimage import map_coordinates

        inverse = np.linalg.inv(self.transform)
        indices = np.asarray(points) @ inverse[:3, :3].T + inverse[:3, 3]
        return map_coordinates(
            self.values, indices.T, order=1, mode="constant", cval=0.0
        )

    def level(self, value):
        if isinstance(value, dict):
            if value.get("kind") == "relative":
                return float(
                    self.values.mean()
                    + float(value.get("relativeValue", 1)) * self.values.std()
                )
            return float(value.get("absoluteValue", 1))
        return float(value)


def read_grid(path, params=None):
    params = params or {}
    path = local_path(path)
    suffix = path.suffix.lower()
    if suffix in (".ccp4", ".map", ".mrc"):
        return read_ccp4(path.read_bytes(), path.name)
    if suffix in (".cube", ".cub"):
        rows = path.read_text().splitlines()
        header = rows[2].split()
        count = int(header[0])
        origin = np.array(header[1:4], float)
        axes = [rows[i].split() for i in range(3, 6)]
        dims = [abs(int(a[0])) for a in axes]
        angstrom = all(int(a[0]) < 0 for a in axes)
        transform = np.eye(4)
        transform[:3, :3] = np.array([a[1:4] for a in axes], float).T
        transform[:3, 3] = origin
        if not angstrom:
            transform[:3] *= 0.529177210859
        start = 6 + abs(count)
        tokens = " ".join(rows[start:]).split()
        orbitals = 1
        if count < 0:
            orbitals = int(tokens[0])
            tokens = tokens[orbitals + 1 :]
        values = np.array(tokens, float).reshape(*dims, orbitals)
        index = int(params.get("orbitalIndex", 0))
        if not 0 <= index < orbitals:
            raise ValueError("orbitalIndex is outside the Cube dataset")
        return Grid(values[..., index], transform, path.name)
    if suffix == ".dx":
        text = path.read_text()
        dims = re.search(r"gridpositions counts\s+(\d+)\s+(\d+)\s+(\d+)", text)
        origin = re.search(r"^origin\s+(.+)$", text, re.MULTILINE)
        axes = re.findall(r"^delta\s+(.+)$", text, re.MULTILINE)
        numbers = re.search(
            r"data follows\s*([\s\S]+?)(?:\nattribute|\nobject|\ncomponent|$)", text
        )
        if not dims or not origin or len(axes) != 3 or not numbers:
            raise ValueError(
                "DX requires counts, origin, three deltas, and scalar data"
            )
        shape = tuple(map(int, dims.groups()))
        values = np.fromstring(numbers[1], sep=" ")
        transform = np.eye(4)
        transform[:3, :3] = np.array([a.split() for a in axes], float).T
        transform[:3, 3] = np.array(origin[1].split(), float)
        return Grid(values[: np.prod(shape)].reshape(shape), transform, path.name)
    if suffix == ".npz":
        with np.load(path, allow_pickle=False) as archive:
            values = archive["values"].copy()
            transform = (
                archive["transform"].copy() if "transform" in archive else np.eye(4)
            )
        return Grid(values, transform, path.name)
    if suffix in (".cif", ".bcif"):
        tables = read_cif(path)
        info = tables.get("volume_data_3d_info", [])
        data = tables.get("volume_data_3d", [])
        if not info or not data:
            raise ValueError(
                "CIF volume requires volume_data_3d_info and volume_data_3d"
            )
        import gemmi

        row = info[0]
        cell = gemmi.UnitCell(
            *(float(row[f"spacegroup_cell_size[{i}]"]) for i in range(3)),
            *(float(row[f"spacegroup_cell_angles[{i}]"]) for i in range(3)),
        )
        count = np.array([int(row[f"sample_count[{i}]"]) for i in range(3)])
        order = [int(row[f"axis_order[{i}]"]) for i in range(3)]
        scale = np.array([float(row[f"dimensions[{i}]"]) for i in range(3)]) / count
        origin = np.array([float(row[f"origin[{i}]"]) for i in range(3)])
        transform = np.eye(4)
        transform[:3, :3] = np.asarray(cell.orth.mat) @ np.diag(scale)
        transform[:3, 3] = np.asarray(cell.orth.mat) @ origin
        values = (
            np.array([r["values"] for r in data], np.float32)
            .reshape(tuple(count[::-1]))
            .transpose(2, 1, 0)
        )
        if order != [0, 1, 2]:
            values = values.transpose(np.argsort(order))
        return Grid(values, transform, path.name, {"tables": tables})
    raise ValueError(f"Unsupported grid file type: {suffix}")


def read_ccp4(raw, label="map"):
    """Decode CCP4/MRC without normalizing values or discarding a nonorthogonal cell."""
    if len(raw) < 1024:
        raise ValueError("CCP4/MRC header is truncated")
    endian = "<"
    header = np.frombuffer(raw[:1024], dtype="<i4")
    if any(n < 1 or n > 1_000_000 for n in header[:3]):
        endian = ">"
        header = np.frombuffer(raw[:1024], dtype=">i4")
    floats = np.frombuffer(raw[:1024], dtype=endian + "f4")
    nx, ny, nz, mode = map(int, header[:4])
    types = {0: "i1", 1: "i2", 2: "f4", 6: "u2", 12: "f2"}
    if mode not in types or min(nx, ny, nz) <= 0:
        raise ValueError(f"Unsupported CCP4/MRC mode {mode}")
    values = (
        np.frombuffer(
            raw,
            dtype=endian + types[mode],
            count=nx * ny * nz,
            offset=1024 + int(header[23]),
        )
        .reshape(nz, ny, nx)
        .transpose(2, 1, 0)
    )
    axes = np.array(header[16:19], int) - 1
    if set(axes) != {0, 1, 2}:
        raise ValueError("CCP4/MRC MAPC/MAPR/MAPS must be a permutation of 1,2,3")
    import gemmi

    lengths = floats[10:13]
    angles = floats[13:16]
    cell = gemmi.UnitCell(*lengths, *(angles if np.all(angles > 0) else (90, 90, 90)))
    divisions = np.maximum(header[7:10], 1)
    basis = np.asarray(cell.orth.mat) / divisions[None, :]
    transform = np.eye(4)
    transform[:3, :3] = basis[:, axes]
    transform[:3, 3] = transform[:3, :3] @ header[4:7]
    if np.any(floats[49:52]):
        transform[:3, 3] = floats[49:52]
    return Grid(values.copy(), transform, label, {"spacegroup": int(header[22])})


def pymol_grid(cmd, name, state=1):
    """Read the PyMOL 3.1 map session layout, including its index-to-world transform."""
    session = cmd.get_session(name, partial=1)
    entries = [
        row
        for row in session.get("names", [])
        if row and row[0] == name and row[4] == 2
    ]
    if not entries:
        raise ValueError(f"{name!r} is not a PyMOL map object")
    states = entries[0][5][2]
    if not 1 <= state <= len(states) or states[state - 1] is None:
        raise ValueError(f"Map {name!r} has no state {state}")
    data = states[state - 1]
    shape = tuple(map(int, data[14][0]))
    count = int(np.prod(shape))

    def values_in(node):
        if isinstance(node, (bytes, bytearray)) and len(node) == count * 4:
            return np.frombuffer(node, dtype=np.float32).copy()
        if isinstance(node, (list, tuple)):
            if len(node) == count and all(isinstance(x, (float, int)) for x in node):
                return np.array(node, np.float32)
            for child in node:
                value = values_in(child)
                if value is not None:
                    return value
        return None

    values = values_in(data[14][2])
    if values is None:
        raise ValueError(
            "Unsupported PyMOL map field layout; provide the original local map file"
        )
    transform = np.eye(4)
    if data[5] is not None:
        transform[:3, :3] = np.diag(data[5])
        transform[:3, 3] = data[2]
    else:
        import gemmi

        lengths, angles = data[1][0]
        cell = gemmi.UnitCell(*lengths, *angles)
        transform[:3, :3] = np.asarray(cell.orth.mat) / np.asarray(data[10])[None, :]
        transform[:3, 3] = transform[:3, :3] @ np.array(data[11])
    matrix = cmd.get_object_matrix(name, state=state)
    if matrix is not None:
        transform = np.array(matrix).reshape(4, 4) @ transform
    return Grid(values.reshape(shape), transform, name)


def _decode_bcif(data, encodings):
    for encoding in reversed(encodings):
        kind = encoding["kind"]
        if kind == "ByteArray":
            dtype = {
                1: "i1",
                2: "<i2",
                3: "<i4",
                4: "u1",
                5: "<u2",
                6: "<u4",
                32: "<f4",
                33: "<f8",
            }[encoding["type"]]
            data = np.frombuffer(data, dtype=dtype).copy()
        elif kind == "FixedPoint":
            data = np.asarray(data, float) / encoding["factor"]
        elif kind == "IntervalQuantization":
            data = encoding["min"] + np.asarray(data) * (
                (encoding["max"] - encoding["min"]) / (encoding["numSteps"] - 1)
            )
        elif kind == "RunLength":
            values = np.asarray(data).reshape(-1, 2)
            data = np.repeat(values[:, 0], values[:, 1].astype(int))
        elif kind == "Delta":
            data = np.cumsum(data, dtype=np.int64) + encoding["origin"]
        elif kind == "IntegerPacking":
            limit = (
                2 ** (8 * encoding["byteCount"] - (0 if encoding["isUnsigned"] else 1))
                - 1
            )
            low = 0 if encoding["isUnsigned"] else -limit - 1
            out, value = [], 0
            for item in data:
                value += int(item)
                if item != limit and (encoding["isUnsigned"] or item != low):
                    out.append(value)
                    value = 0
            data = np.array(out)
        elif kind == "StringArray":
            offsets = _decode_bcif(encoding["offsets"], encoding["offsetEncoding"])
            indices = _decode_bcif(data, encoding["dataEncoding"])
            strings = [
                encoding["stringData"][a:b] for a, b in itertools.pairwise(offsets)
            ]
            data = [strings[int(i)] if i >= 0 else "" for i in indices]
        else:
            raise ValueError(f"Unsupported BinaryCIF encoding: {kind}")
    return data


def read_cif(path):
    path = local_path(path)
    if path.suffix.lower() == ".bcif":
        import msgpack

        blob = msgpack.unpackb(path.read_bytes(), raw=False)
        result = {}
        for block in blob["dataBlocks"]:
            for category in block["categories"]:
                columns = {}
                for col in category["columns"]:
                    values = list(
                        _decode_bcif(col["data"]["data"], col["data"]["encoding"])
                    )
                    if col.get("mask"):
                        mask = _decode_bcif(
                            col["mask"]["data"], col["mask"]["encoding"]
                        )
                        values = [v if m == 0 else None for v, m in zip(values, mask)]
                    columns[col["name"]] = values
                result[category["name"].lstrip("_")] = [
                    dict(zip(columns, row)) for row in zip(*columns.values())
                ]
        return result
    import gemmi

    text = (
        gzip.decompress(path.read_bytes()).decode()
        if path.suffix == ".gz"
        else path.read_text()
    )
    doc = gemmi.cif.read_string(text)
    result = {}
    for block in doc:
        for prefix in block.get_mmcif_category_names():
            cols = block.get_mmcif_category(prefix)
            result[prefix.strip("_.")] = [
                dict(zip(cols, row)) for row in zip(*cols.values())
            ]
    return result


def read_data(value):
    if value is None or (isinstance(value, str) and not value):
        return {}
    if isinstance(value, dict):
        return value.copy()
    path = local_path(value)
    suffix = path.suffix.lower()
    if suffix in (".map", ".ccp4", ".mrc", ".cube", ".cub", ".dx"):
        return {"grid": read_grid(path), "file": str(path)}
    if suffix == ".npz":
        with np.load(path, allow_pickle=False) as archive:
            content = {k: archive[k].copy() for k in archive}
        if "values" in content and np.asarray(content["values"]).ndim == 3:
            content["grid"] = Grid(
                content.pop("values"), content.pop("transform", np.eye(4)), path.name
            )
        return content
    if suffix in (".cif", ".bcif") or path.name.endswith(".cif.gz"):
        return {"tables": read_cif(path), "file": str(path)}
    if suffix in (".json", ".mvsj"):
        content = json.loads(path.read_text())
        if not isinstance(content, dict):
            if (
                isinstance(content, list)
                and len(content) == 1
                and isinstance(content[0], dict)
            ):
                content = content[0]
            else:
                raise ValueError("Data JSON must contain an object")
        content = {**content, "_base": str(path.parent)}
        if isinstance(content.get("grid"), str):
            content["grid"] = read_grid(local_path(content["grid"], path.parent))
        elif isinstance(content.get("grid"), dict):
            g = content["grid"]
            content["grid"] = Grid(
                g["values"], g.get("transform", np.eye(4)), g.get("label", "grid")
            )
        return content
    return {"file": str(path), "_base": str(path.parent)}


def atom_fields(atoms, data):
    """Align explicit annotations by atom/residue identity, never by partial row order."""
    result = {}
    count = len(atoms)
    for key, values in data.get("atom_data", {}).items():
        if len(values) != count:
            raise ValueError(
                f"atom_data.{key} has {len(values)} entries; selection contains {count} atoms"
            )
        result[key] = np.asarray(values)
    rows = data.get("residues", [])
    for row in rows:
        if "resi" not in row or "chain" not in row:
            raise ValueError("Residue annotations require chain and resi")
    if rows:
        properties = set().union(*(r.keys() for r in rows)) - {
            "chain",
            "resi",
            "model",
            "segi",
        }
        for key in properties:
            values = []
            for atom in atoms:
                matches = [
                    r
                    for r in rows
                    if str(r["resi"]) == atom.resi
                    and str(r["chain"]) == atom.chain
                    and r.get("model", atom.model) == atom.model
                    and r.get("segi", atom.segi) == atom.segi
                ]
                if len(matches) > 1:
                    raise ValueError(
                        f"Ambiguous annotation for chain {atom.chain} residue {atom.resi}"
                    )
                values.append(matches[0].get(key) if matches else None)
            result[key] = np.array(values, object)
    return result
