"""Exercise real Qt/OpenGL and native ray renderers in an isolated process."""

import argparse
import json
from pathlib import Path
from time import monotonic, perf_counter, sleep

import numpy as np
from PIL import Image
from samples import data_for, molecule


def composited_pixels(path):
    image = Image.open(path).convert("RGBA")
    base = Image.new("RGBA", image.size, "white")
    base.alpha_composite(image)
    return np.asarray(base.convert("RGB"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gui", action="store_true")
    parser.add_argument("--only", nargs="+")
    parser.add_argument("--benchmark", action="store_true")
    parser.add_argument("--structure", type=Path)
    parser.add_argument("--visuals", action="store_true")
    parser.add_argument("--output", type=Path, default=Path(".cache/validation"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.gui:
        import pymol

        pymol.invocation.options.show_splash = 0
        from pmg_qt.pymol_qt_gui import PyMOLApplication, PyMOLQtGUI

        app = PyMOLApplication(["molstar-validation"])
        window = PyMOLQtGUI()
        widget = window.pymolwidget
        cmd = widget.cmd

        def in_context(function):
            with widget:
                return function()

        cmd._call_with_opengl_context = in_context

        def pump(seconds=0.08):
            end = monotonic() + seconds
            while monotonic() < end:
                app.processEvents()
                sleep(0.002)

        window.resize(800, 650)
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
    else:
        import pymol2

        instance = pymol2.PyMOL()
        instance.start()
        cmd = instance.cmd

        def pump(seconds=0):
            pass

    from molstar_style_in_pymol import molstar_style
    from molstar_style_in_pymol.controller import manager_for
    from molstar_style_in_pymol.registry import (
        EFFECT_STYLES,
        EXTENSION,
        MATERIALS,
        PARTICLE,
        SHAPE,
        STRUCTURE,
        VOLUME,
    )

    state = molecule()
    cmd.load_model(state.model, "sample")
    if args.structure:
        cmd.delete("sample")
        cmd.load(str(args.structure), "sample")
        cmd.dss("sample")
    cmd.hide("everything", "sample")
    cmd.orient("sample")
    cmd.bg_color("white")
    cmd.set("orthoscopic", 1)
    original = cmd.get_coords("sample").copy()
    styles = args.only or [
        *STRUCTURE,
        *VOLUME,
        *PARTICLE,
        *SHAPE,
        *EXTENSION,
        *MATERIALS,
        *EFFECT_STYLES,
    ]
    if args.benchmark:
        styles = []
    if args.visuals:
        from molstar_style_in_pymol.registry import reference

        styles = [
            f"{rep}:{visual}"
            for rep, visuals in reference()["visuals"].items()
            for visual in visuals
        ]
    report = {}
    failed = {}
    for rep in styles:
        label = rep.replace(":", "__")
        rep, _, visual = rep.partition(":")
        print("Rendering", rep, flush=True)
        t = perf_counter()
        try:
            data = data_for(rep, state)
            if rep == "mvs":
                model_file = args.output / "mvs_model.pdb"
                cmd.save(str(model_file), "sample")
                data = {
                    "_base": str(args.output),
                    "root": {
                        "kind": "root",
                        "children": [
                            {
                                "kind": "download",
                                "params": {"url": model_file.name},
                                "children": [
                                    {
                                        "kind": "parse",
                                        "params": {"format": "pdb"},
                                        "children": [
                                            {
                                                "kind": "structure",
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
            entry = molstar_style(
                rep,
                "sample",
                quality="low",
                data=data,
                params={"resolution": 1.2, **({"visuals": [visual]} if visual else {})},
                quiet=1,
                _self=cmd,
            )
            cmd.orient(entry.name)
            cmd.zoom(entry.name, 2)
            cmd.clip("slab", 200)
            pump()
            assert entry.name in manager_for(cmd).entries, "OpenGL callback failed"
            if args.gui:
                filename = args.output / f"{label}-gpu.png"
                molstar_style(
                    "png", filename=str(filename), width=480, height=360, _self=cmd
                )
                pixels = composited_pixels(filename)
                assert np.count_nonzero(pixels.min(axis=2) < 235) > 30, (
                    "Empty GPU image"
                )
                errors = [
                    d.error for ds in entry.drawings.values() for d in ds if d.error
                ]
                assert not errors, errors
            filename = args.output / f"{label}-ray.png"
            if rep == "pairwise-metric":
                from molstar_style_in_pymol.panels import raster

                raster(entry.panels[0], 360).save(filename)
            else:
                molstar_style(
                    "ray", filename=str(filename), width=240, height=180, _self=cmd
                )
            pixels = composited_pixels(filename)
            assert np.count_nonzero(pixels.min(axis=2) < 235) > 10, "Empty ray image"
            molstar_style("reset", _self=cmd)
            np.testing.assert_array_equal(cmd.get_coords("sample"), original)
            report[label] = {
                "seconds": round(perf_counter() - t, 3),
                "geometry_bytes": entry.nbytes,
                "gpu": args.gui,
                "ray": True,
            }
        except Exception as exc:
            import traceback

            traceback.print_exc()
            failed[label] = str(exc)
            manager_for(cmd).reset("all")
    if args.benchmark:
        if not args.gui:
            raise ValueError("Benchmark requires --gui")
        from benchmark import benchmark

        benchmark(cmd, widget, pump, args.output, molstar_style, report, args.structure)
    (args.output / "report.json").write_text(
        json.dumps({"passed": report, "failed": failed}, indent=2) + "\n"
    )
    if args.gui:
        window.close()
        pump()
    else:
        instance.stop()
    print(json.dumps({"passed": len(report), "failed": failed}, indent=2), flush=True)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
