"""Render PyMOL at the exact camera exported by the Mol* comparison viewer."""

import argparse
import json
from pathlib import Path
from time import monotonic, sleep

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path(".cache/reference"))
    parser.add_argument("--name", default="pymol")
    parser.add_argument(
        "--source",
        type=Path,
        help="Optional isolated package source for a before render",
    )
    parser.add_argument("--native-ss", action="store_true")
    parser.add_argument("--occlusion", action="store_true")
    args = parser.parse_args()
    if args.source:
        import sys

        sys.path.insert(0, str(args.source))
    import pymol

    from molstar_style_in_pymol import molstar_style
    from molstar_style_in_pymol.controller import manager_for

    pymol.invocation.options.show_splash = 0
    from pmg_qt.pymol_qt_gui import PyMOLApplication, PyMOLQtGUI

    app = PyMOLApplication(["molstar-reference-comparison"])
    window = PyMOLQtGUI()
    widget, cmd = window.pymolwidget, window.pymolwidget.cmd

    def context(function):
        with widget:
            return function()

    def pump(seconds=0.2):
        deadline = monotonic() + seconds
        while monotonic() < deadline:
            app.processEvents()
            sleep(0.002)

    cmd._call_with_opengl_context = context
    window.show()
    for key in (
        "internal_gui",
        "internal_feedback",
        "internal_prompt",
        "seq_view",
        "movie_panel",
    ):
        cmd.set(key, 0)
    ratio = widget.devicePixelRatioF()
    widget.setFixedSize(round(1200 / ratio), round(900 / ratio))
    window.adjustSize()
    pump()
    with widget:
        widget.resizeGL(widget.width(), widget.height())
    ref = json.loads(args.reference.read_text())
    cmd.load(str(args.structure), "sample")
    keys = {(s["chain"], s["resi"]) for s in ref["segments"]}
    indices = [
        a.index for a in cmd.get_model("sample").atom if (a.chain, a.resi) in keys
    ]
    cmd.select_list("reference_atoms", "sample", indices, mode="index")
    if not args.native_ss:
        ss = {(s["chain"], s["resi"]): s["ss"] for s in ref["segments"]}
        cmd.alter(
            "sample", "ss = assignment.get((chain, resi), '')", space={"assignment": ss}
        )
    cmd.hide("everything")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("antialias", 2)
    camera = ref["camera"]
    matrix = np.asarray(ref["view"]).reshape(4, 4).T
    distance = np.linalg.norm(np.array(camera["position"]) - camera["target"])
    cmd.set("field_of_view", np.rad2deg(camera["fov"]))
    cmd.set("orthoscopic", int(camera["mode"] == "orthographic"))
    cmd.set_view(
        [
            *matrix[:3, :3].T.ravel(),
            0,
            0,
            -distance,
            *camera["target"],
            max(1, distance - camera["radiusMax"]),
            distance + camera["radiusMax"],
            np.rad2deg(camera["fov"]),
        ]
    )
    cmd.set("orthoscopic", int(camera["mode"] == "orthographic"))
    options = ref.get("options", {})
    params = {
        "linearSegments": 8,
        "radialSegments": 16,
        "visuals": ["polymer-trace"],
        "postprocessing": {"occlusion": args.occlusion},
        "colorParams": {"value": 0x1B9E77},
        **options.get("params", {}),
    }
    if "value" not in params["colorParams"]:
        params["colorParams"]["value"] = 0x1B9E77
    args.output.mkdir(parents=True, exist_ok=True)
    try:
        molstar_style(
            "cartoon",
            "reference_atoms and polymer",
            color=options.get("color", "uniform"),
            params=params,
            _self=cmd,
        )
        pump()
        drawing = next(iter(manager_for(cmd).active_drawings()))
        for renderer, operation in (("gpu", "png"), ("ray", "ray")):
            molstar_style(
                operation,
                filename=str(args.output / f"{args.name}-{renderer}.png"),
                width=1200,
                height=900,
                _self=cmd,
            )
        modelview, projection, viewport = drawing.matrices
        report = {
            "view": list(cmd.get_view()),
            "modelview": modelview.tolist(),
            "projection": projection.tolist(),
            "viewport": viewport.tolist(),
            "error": drawing.error,
        }
        (args.output / f"{args.name}-camera.json").write_text(json.dumps(report) + "\n")
    finally:
        molstar_style("reset", name="all", quiet=1, _self=cmd)
        cmd.delete("all")
        window.hide()
        pump()


if __name__ == "__main__":
    main()
