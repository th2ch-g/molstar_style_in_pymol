"""Render the documentation gallery in a separate real PyMOL Qt process."""

import argparse
import hashlib
import json
import subprocess
from pathlib import Path
from time import monotonic, sleep

import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]

from molstar_style_in_pymol import molstar_style  # noqa: E402
from molstar_style_in_pymol.controller import manager_for  # noqa: E402
from molstar_style_in_pymol.registry import (  # noqa: E402
    EFFECT_STYLES,
    EXTENSION,
    MATERIALS,
    PARTICLE,
    PRESETS,
    SHAPE,
    STRUCTURE,
    VOLUME,
)
from molstar_style_in_pymol.source import read  # noqa: E402

GROUPS = {
    "Structure": (*STRUCTURE, "nucleic"),
    "Volume": VOLUME,
    "Particles": PARTICLE,
    "Measurements and shapes": SHAPE,
    "Extensions": EXTENSION,
    "Materials": tuple(MATERIALS),
    "Effects": EFFECT_STYLES,
    "Presets": tuple(p for p in PRESETS if p != "empty"),
}
ATOMIC = {"ball-and-stick", "ellipsoid", "line", "point", "spacefill", "label"}
SYNTHETIC = {*VOLUME, *PARTICLE, *SHAPE, *EXTENSION, "carbohydrate", "polyhedron"}


def specimen(style):
    if style == "nucleic":
        return "DNA built with PyMOL fnab"
    if style in ATOMIC:
        return "PyMOL tryptophan fragment"
    if style in SYNTHETIC and style != "mvs":
        return "Synthetic demonstration data"
    if style.startswith(("quality-", "validation-")) or style == "partial-charges":
        return "1CRN with synthetic annotation values"
    return "Crambin (PDB 1CRN)"


def mvs_input(structure):
    representation = {"kind": "representation", "params": {"type": "cartoon"}}
    component = {
        "kind": "component",
        "params": {"selector": "all"},
        "children": [representation],
    }
    model = {"kind": "structure", "children": [component]}
    parse = {
        "kind": "parse",
        "params": {"format": "mmcif" if structure.suffix.lower() == ".cif" else "pdb"},
        "children": [model],
    }
    download = {
        "kind": "download",
        "params": {"url": structure.name},
        "children": [parse],
    }
    return {
        "_base": str(structure.parent),
        "root": {"kind": "root", "children": [download]},
    }


def load_sample(cmd, style, structure):
    # Initialize PyMOL before importing chempy so its fragment paths are set.
    from samples import data_for, molecule

    cmd.delete("all")
    if style == "nucleic":
        cmd.fnab("ATGCGCAT", name="sample")
    elif style in ATOMIC:
        cmd.fragment("trp", "sample")
        if style == "ellipsoid":
            model = cmd.get_model("sample")
            for atom in model.atom:
                atom.u_aniso = [0.12, 0.3, 0.6, 0.03, 0.02, 0.02]
            cmd.delete("sample")
            cmd.load_model(model, "sample")
    elif style in SYNTHETIC and style != "mvs":
        cmd.load_model(molecule().model, "sample")
        if style == "carbohydrate":
            cmd.remove("sample and not resn NAG+MAN")
        elif style == "polyhedron":
            cmd.remove("sample and not chain Z+W")
    else:
        cmd.load(str(structure), "sample")
        cmd.remove("sample and not polymer.protein")
        cmd.dss("sample")
    cmd.hide("everything", "all")
    cmd.bg_color("white")
    cmd.orient("sample")
    cmd.turn("y", 20)
    cmd.turn("z", -15)
    cmd.zoom("sample", 3, complete=1)
    cmd.clip("slab", 200)
    snapshots, _ = read(cmd, "sample")
    state = next(iter(snapshots.values()))[0]
    return mvs_input(structure) if style == "mvs" else data_for(style, state)


def validate_image(path, width, height):
    with Image.open(path) as image:
        image.load()
        rgba = image.convert("RGBA")
    white = Image.new("RGBA", rgba.size, "white")
    white.alpha_composite(rgba)
    pixels = np.asarray(white.convert("RGB"))
    visible = int(np.count_nonzero(pixels.min(axis=2) < 235))
    if visible < 100:
        raise ValueError(f"Empty gallery image: {path.name}")
    # Metric panels are square; center them on the same documentation canvas.
    image = white.convert("RGB")
    if image.size != (width, height):
        image.thumbnail((width, height), Image.Resampling.LANCZOS)
        canvas = Image.new("RGB", (width, height), "white")
        canvas.paste(image, ((width - image.width) // 2, (height - image.height) // 2))
        image = canvas
    image.save(path, optimize=True)
    return {
        "visible_pixels": visible,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--structure", type=Path, required=True, help="Local 1CRN CIF or PDB"
    )
    parser.add_argument("--output", type=Path, default=Path("docs/gallery"))
    parser.add_argument(
        "--manifest", type=Path, default=Path(".cache/gallery/manifest.json")
    )
    parser.add_argument("--only", nargs="+")
    parser.add_argument("--width", type=int, default=1200)
    parser.add_argument("--height", type=int, default=900)
    parser.add_argument("--quality", default="high")
    args = parser.parse_args()
    styles = [style for group in GROUPS.values() for style in group]
    if args.only:
        unknown = set(args.only) - set(styles)
        if unknown:
            parser.error(f"Unknown gallery entries: {sorted(unknown)}")
        styles = [style for style in styles if style in args.only]
    structure = args.structure.resolve()
    if not structure.is_file():
        parser.error("--structure must be an existing local CIF or PDB file")
    if structure.suffix.lower() not in (".cif", ".pdb"):
        parser.error("--structure must have a .cif or .pdb suffix")
    if min(args.width, args.height) < 100:
        parser.error("Image dimensions must be at least 100 pixels")
    args.output.mkdir(parents=True, exist_ok=True)

    import pymol

    pymol.invocation.options.show_splash = 0
    from pmg_qt.pymol_qt_gui import PyMOLApplication, PyMOLQtGUI

    app = PyMOLApplication(["molstar-gallery"])
    window = PyMOLQtGUI()
    widget, cmd = window.pymolwidget, window.pymolwidget.cmd

    def in_context(function):
        with widget:
            return function()

    def pump(seconds=0.08):
        end = monotonic() + seconds
        while monotonic() < end:
            app.processEvents()
            sleep(0.002)

    cmd._call_with_opengl_context = in_context
    window.resize(1000, 800)
    window.show()
    pump(0.3)
    for key in (
        "internal_gui",
        "internal_feedback",
        "internal_prompt",
        "seq_view",
        "movie_panel",
    ):
        cmd.set(key, 0)
    ratio = widget.devicePixelRatioF()
    widget.setFixedSize(round(args.width / ratio), round(args.height / ratio))
    window.adjustSize()
    pump(0.3)
    with widget:
        widget.resizeGL(widget.width(), widget.height())
    cmd.set("orthoscopic", 1)
    cmd.set("antialias", 2)
    cmd.set("ray_opaque_background", 1)
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    report = {
        "revision": revision,
        "structure_sha256": hashlib.sha256(structure.read_bytes()).hexdigest(),
        "quality": args.quality,
        "width": args.width,
        "height": args.height,
        "images": {},
    }
    try:
        for style in styles:
            print(f"Rendering {style}", flush=True)
            data = load_sample(cmd, style, structure)
            entry = molstar_style(
                style,
                "sample",
                quality=args.quality,
                data=data,
                quiet=1,
                _self=cmd,
            )
            if style in SYNTHETIC and style not in ("mvs", "pairwise-metric"):
                cmd.orient(entry.name)
                cmd.turn("y", 20)
                cmd.turn("x", 15)
                cmd.zoom(entry.name, 2, complete=1)
                if style == "direct-volume":
                    view = list(cmd.get_view())
                    view[11] *= 0.45
                    cmd.set_view(view)
                elif style in ("orbital", "orbital-density"):
                    cmd.turn("y", 50)
                cmd.clip("slab", 200)
            pump()
            if entry.name not in manager_for(cmd).entries:
                raise RuntimeError(f"OpenGL callback failed for {style}")
            record = {"specimen": specimen(style), "command": f"molstar_style {style}"}
            for renderer, operation in (("gpu", "png"), ("ray", "ray")):
                path = args.output / f"{style}-{renderer}.png"
                molstar_style(
                    operation,
                    filename=str(path),
                    width=args.width,
                    height=args.height,
                    _self=cmd,
                )
                record[renderer] = validate_image(path, args.width, args.height)
            report["images"][style] = record
            molstar_style("reset", name="all", quiet=1, _self=cmd)
    finally:
        molstar_style("reset", name="all", quiet=1, _self=cmd)
        cmd.delete("all")
        window.hide()
        pump()
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(report, indent=2) + "\n")
    print(
        f"Validated {len(report['images'])} GPU/ray pairs in {args.output}", flush=True
    )


if __name__ == "__main__":
    main()
