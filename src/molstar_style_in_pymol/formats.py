"""Local extension formats and model annotation normalization."""

import re
import zlib

import numpy as np

from .data import finite, local_path


def hydrate(data, atoms):
    data = dict(data)
    tables = data.get("tables", {})
    if not tables:
        return data
    atom_rows = tables.get("atom_site", [])
    byres = {}
    label_to_auth = {}
    for row in atom_rows:
        chain = str(row.get("auth_asym_id", row.get("label_asym_id", "")))
        resi = str(row.get("auth_seq_id", row.get("label_seq_id", "")))
        ins = row.get("pdbx_PDB_ins_code")
        if ins and ins not in (".", "?"):
            resi += str(ins)
        label_to_auth[(str(row.get("label_asym_id")), str(row.get("label_seq_id")))] = (
            chain,
            resi,
        )
        byres.setdefault(
            (chain, resi),
            {"chain": chain, "resi": resi, "entity_id": row.get("label_entity_id")},
        )
    metrics = {
        str(r["id"]): str(r.get("name", r.get("type", ""))).lower()
        for r in tables.get("ma_qa_metric", [])
    }
    for row in tables.get("ma_qa_metric_local", []):
        key = label_to_auth.get((str(row["label_asym_id"]), str(row["label_seq_id"])))
        if key:
            name = metrics.get(str(row["metric_id"]), "")
            field = (
                "plddt"
                if "plddt" in name
                else "qmean"
                if "qmean" in name
                else name.replace(" ", "_")
            )
            byres[key][field] = float(row["metric_value"])
    sources = {}
    for name in ("entity_src_gen", "entity_src_nat", "pdbx_entity_src_syn"):
        for row in tables.get(name, []):
            fields = [
                v
                for k, v in row.items()
                if "scientific_name" in k or "taxonomy_id" in k
            ]
            sources[str(row["entity_id"])] = " / ".join(str(v) for v in fields)
    for row in byres.values():
        if str(row["entity_id"]) in sources:
            row["entity_source"] = sources[str(row["entity_id"])]
    data["residues"] = [*byres.values(), *data.get("residues", [])]
    for category, field, values in (
        ("pdbx_validate_rmsd_bond", "geometry_quality", ("Z",)),
        ("pdbx_validate_rmsd_angle", "geometry_quality", ("Z",)),
        ("pdbx_vrpt_model_instance_map_fitting", "density_fit", ("RSCC",)),
        ("pdbx_nmr_software", "random_coil_index", ("rci",)),
    ):
        for row in tables.get(category, []):
            chain = str(row.get("auth_asym_id", row.get("auth_asym_id_1", "")))
            resi = str(row.get("auth_seq_id", row.get("auth_seq_id_1", "")))
            value = next(
                (row[k] for k in values if row.get(k) not in (None, ".", "?")), None
            )
            if value is not None and (chain, resi) in byres:
                byres[(chain, resi)][field] = abs(float(value))
    if "cell" in tables:
        row = tables["cell"][0]
        data["cell"] = [
            float(row[k])
            for k in (
                "length_a",
                "length_b",
                "length_c",
                "angle_alpha",
                "angle_beta",
                "angle_gamma",
            )
        ]
    return data


def g3d(data, p):
    import msgpack

    if "file" in data:
        path = local_path(data["file"])
        with path.open("rb") as stream:
            header = msgpack.unpackb(stream.read(64000).rstrip(b"\0"), raw=False)
            if header.get("magic") != "G3D":
                raise ValueError("Invalid G3D magic")
            resolution = int(p.get("resolution", min(header["resolutions"])))
            block = header["offsets"][str(resolution)]
            stream.seek(block["offset"])
            payload = stream.read(block["size"])
        try:
            payload = zlib.decompress(payload)
        except zlib.error:
            payload = zlib.decompress(payload, -15)
        groups = msgpack.unpackb(payload, raw=False)
    else:
        groups = data.get("data", data.get("haplotypes", {}))
    particles = []
    for haplotype, chromosomes in groups.items():
        for chromosome, row in chromosomes.items():
            points = finite(
                np.c_[row["x"], row["y"], row["z"]], (None, 3), "G3D coordinates"
            )
            starts = np.asarray(row["start"])
            mask = np.ones(len(points), bool)
            if p.get("chromosome"):
                mask &= chromosome == p["chromosome"]
            if p.get("haplotype"):
                mask &= haplotype == p["haplotype"]
            if "region" in p:
                mask &= (starts >= p["region"][0]) & (starts <= p["region"][1])
            if mask.any():
                particles.append(
                    {
                        "points": points[mask],
                        "position": points[mask][0],
                        "radius": float(p.get("radius", 0.3)),
                        "entity": chromosome,
                        "label": f"{haplotype} {chromosome}",
                    }
                )
    if not particles:
        raise ValueError("No G3D coordinates match the requested region")
    from .shapes import particles as geometry

    return geometry({"particles": particles}, "particle-fibers", p)


def kinemage(data, p):
    from .primitives import cylinder, indexed, sphere, text_mesh
    from .scene import Geometry, PickTarget
    from .themes import rgb

    source = data.get("text") or local_path(data["file"]).read_text()
    result = Geometry()
    kind = None
    color = rgb(0xFFFFFF)
    previous = None
    triangle = []
    color_alias = {
        "sea": "seagreen",
        "sky": "skyblue",
        "hotpink": "hotpink",
        "lilac": "orchid",
        "peach": "peachpuff",
        "deadwhite": "white",
        "deadblack": "black",
        "invisible": "black",
        "greentint": "palegreen",
        "pinktint": "pink",
        "bluetint": "lightblue",
        "gray": "gray",
    }
    for line in source.splitlines():
        line = line.strip()
        if line.startswith("@"):
            directive = line.split()[0][1:]
            if directive.endswith("list"):
                kind = directive[:-4]
                previous = None
                triangle = []
                match = re.search(r"color\s*=\s*(\w+)", line)
                if match:
                    color = rgb(color_alias.get(match[1], match[1]))
            continue
        match = re.match(r"\{([^}]+)\}(.*)", line)
        if not match or kind is None:
            continue
        label, body = match.groups()
        tokens = body.replace(",", " ").split()
        numbers = []
        for token in reversed(tokens):
            try:
                numbers.append(float(token))
                if len(numbers) == 3:
                    break
            except ValueError:
                pass
        if len(numbers) != 3:
            raise ValueError(f"Missing Kinemage coordinates for {label}")
        point = np.array(numbers[::-1])
        target = (PickTarget(label=label),)
        if kind in ("vector", "ribbon"):
            if previous is not None and "P" not in tokens:
                result.add(
                    cylinder(
                        previous, point, float(p.get("sizeFactor", 0.06)), color, 0
                    ),
                    target,
                )
            previous = point
        elif kind in ("ball", "sphere", "dot"):
            result.add(
                sphere(
                    point,
                    float(p.get("sizeFactor", 0.2 if kind != "dot" else 0.05)),
                    color,
                    0,
                    8,
                ),
                target,
            )
        elif kind in ("label", "word"):
            result.add(
                text_mesh(label, point, color, float(p.get("textSize", 0.4)), 0),
                target,
                True,
            )
        elif kind == "triangle":
            if "P" in tokens:
                triangle = []
            triangle.append(point)
            if len(triangle) >= 3:
                result.add(indexed(triangle[-3:], [[0, 1, 2]], color), target)
        else:
            raise ValueError(f"Unsupported Kinemage list {kind}")
    return result


def local_structure(path):
    """Create an independent snapshot without loading or modifying PyMOL objects."""
    import gemmi
    from chempy import Atom as ModelAtom
    from chempy import Bond
    from chempy.models import Indexed

    from .source import Atom, State

    structure = gemmi.read_structure(str(local_path(path)))
    model = Indexed()
    atoms = []
    coords = []
    for chain in structure[0]:
        for residue in chain:
            info = gemmi.find_tabulated_residue(residue.name)
            kind = (
                "protein"
                if info.is_amino_acid()
                else "nucleic"
                if info.is_nucleic_acid()
                else "other"
            )
            for a in residue:
                i = len(atoms)
                atoms.append(
                    Atom(
                        "",
                        i + 1,
                        a.name,
                        str(residue.seqid),
                        residue.name,
                        chain.name,
                        "",
                        a.element.name.upper(),
                        "",
                        "" + kind,
                        (0.5, 0.5, 0.5),
                        a.element.vdw_r,
                        str(a.altloc).strip("\0"),
                        a.occ,
                        a.b_iso,
                    )
                )
                coords.append([a.pos.x, a.pos.y, a.pos.z])
                m = ModelAtom()
                m.coord = coords[-1]
                m.name = a.name
                m.index = i + 1
                model.atom.append(m)
    xyz = np.asarray(coords, np.float32)
    from scipy.spatial import cKDTree

    for i, j in cKDTree(xyz).query_pairs(2.3):
        a, b = atoms[i], atoms[j]
        threshold = (
            gemmi.Element(a.element).covalent_r
            + gemmi.Element(b.element).covalent_r
            + 0.3
        )
        if np.linalg.norm(xyz[i] - xyz[j]) < threshold and (
            not a.alt or not b.alt or a.alt == b.alt
        ):
            bond = Bond()
            bond.index = [i, j]
            bond.order = 1
            model.bond.append(bond)
    return State(
        tuple(atoms),
        xyz,
        tuple(tuple(b.index) for b in model.bond),
        model,
        np.zeros(len(atoms), int),
        {},
    )


def mvs(data, p, quality, build):
    from .scene import Geometry

    result = Geometry()
    root = data.get("root")
    if root is None and data.get("kind") == "multiple":
        root = data["snapshots"][int(p.get("snapshot", 0))]["root"]
    if root is None:
        raise ValueError("MVS requires root or local snapshots")
    base = data.get("_base", ".")

    def selector(state, value):
        if value in (None, "all"):
            return state
        if isinstance(value, str):
            mask = [
                a.kind == ("protein" if value == "protein" else "nucleic")
                if value in ("protein", "nucleic")
                else a.kind != "other"
                if value == "polymer"
                else a.resn in ("HOH", "WAT")
                if value == "water"
                else a.kind == "other"
                for a in state.atoms
            ]
        else:
            rows = value if isinstance(value, list) else [value]

            def matches(a, r):
                keys = {
                    "auth_asym_id": a.chain,
                    "label_asym_id": a.chain,
                    "auth_seq_id": a.resi,
                    "label_seq_id": a.resi,
                    "auth_atom_id": a.name,
                    "label_atom_id": a.name,
                    "label_comp_id": a.resn,
                    "type_symbol": a.element,
                }
                unknown = set(r) - set(keys) - {"beg_auth_seq_id", "end_auth_seq_id"}
                if unknown:
                    raise ValueError(
                        f"MVS selector fields require mapped data: {sorted(unknown)}"
                    )
                return all(
                    str(keys[k]) == str(v) for k, v in r.items() if k in keys
                ) and int(r.get("beg_auth_seq_id", -(10**9))) <= int(a.resi) <= int(
                    r.get("end_auth_seq_id", 10**9)
                )

            mask = [any(matches(a, row) for row in rows) for a in state.atoms]
        return state.subset(mask)

    def visit(node, state=None, path=None, transform=None):
        kind = node["kind"]
        options = node.get("params", {})
        children = node.get("children", [])
        if kind == "download":
            path = local_path(options["url"], base)
        elif kind == "parse":
            if options.get("format") not in ("mmcif", "pdb", "bcif"):
                raise ValueError(
                    "MVS structure parse format must be mmcif, pdb, or bcif"
                )
            state = local_structure(path)
        elif kind == "component":
            state = selector(state, options.get("selector", "all"))
        elif kind == "transform":
            rotation = np.array(options.get("rotation", np.eye(3).ravel())).reshape(
                3, 3, order="F"
            )
            state = state.subset(np.ones(len(state.atoms), bool))
            state.coords = state.coords @ rotation.T + options.get(
                "translation", [0, 0, 0]
            )
        elif kind == "representation":
            cp = {}
            color = "auto"
            for child in children:
                if child["kind"] == "color":
                    color = "uniform"
                    cp["value"] = child["params"]["color"]
            localp = {**p, **options, "colorParams": cp}
            opacity = next(
                (
                    float(v["params"]["opacity"])
                    for v in children
                    if v["kind"] == "opacity"
                ),
                1,
            )
            g = build(state, options.get("type", "cartoon"), color, localp, {}, quality)
            for piece in g.pieces:
                piece.mesh.opacity *= opacity
            result.extend(g)
        elif kind == "label":
            from .shapes import measurements

            result.extend(
                measurements(state, "custom-label", p, {"text": options["text"]})
            )
        elif kind not in (
            "root",
            "structure",
            "color",
            "opacity",
            "focus",
            "camera",
            "canvas",
            "tooltip",
        ):
            raise ValueError(
                f"MVS node {kind!r} needs normalized local annotation data"
            )
        if kind in ("camera", "canvas", "focus"):
            result.metadata[kind] = options
        for child in children:
            visit(child, state, path, transform)

    visit(root)
    return result
