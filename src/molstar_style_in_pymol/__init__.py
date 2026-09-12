"""Independent Mol*-inspired visualization; importing does not import PyMOL."""

__version__ = "0.1.0"


def molstar_style(
    style="default",
    selection="all",
    representation=None,
    color="auto",
    quality="medium",
    name="molstar",
    transparency="keep",
    state=0,
    data=None,
    params=None,
    cache_mb=2048,
    gpu_cache_mb=256,
    filename="",
    width=0,
    height=0,
    quiet=0,
    _self=None,
):
    """
    DESCRIPTION
        Apply independent Mol*-inspired geometry and materials in PyMOL 3.1.
        Uses local inputs only. All loaded states are prepared before playback.

    USAGE
        molstar_style [style [, selection [, representation [, color [, quality [, name]]]]]]
        molstar_style list
        molstar_style help
        molstar_style refresh, name=all
        molstar_style reset, name=all
        molstar_style png, filename=figure.png, width=1600, height=1200
        molstar_style ray, filename=figure.png, width=1600, height=1200

    ARGUMENTS
        style: representation, preset, material, or effect; list shows all names.
        selection: PyMOL molecular selection or a local map object.
        representation: explicit geometry override for material/effect styles.
        color: auto, keep, or a Mol* color theme; see list.
        quality: lowest, lower, low, medium, high, higher, highest, auto, custom.
        name: managed group; reapplying replaces it after successful preparation.
        transparency: keep, or a fraction from 0 to 1.
        state: 0 prepares all states; a positive integer selects one state.
        data: local file path or Python dictionary with scientific inputs.
        params: local JSON path or Python dictionary of representation parameters.
        cache_mb, gpu_cache_mb: CPU geometry and GPU buffer budgets in MiB.

    NOTES
        Node.js, Mol*, and CueMol are not runtime dependencies. Native ray uses
        the same meshes; direct volumes use sampled density planes for ray.
        Refresh after editing source coordinates, colors, or annotations.
        Reset restores source representations. See docs/pymol_molstar.md for schemas
        and the per-feature approximation and validation details.
    """
    if _self is None:
        from pymol import cmd as _self
    if not _self.is_gui_thread() and callable(
        getattr(_self, "_call_in_gui_thread", None)
    ):
        arguments = locals().copy()
        return _self._call_in_gui_thread(lambda: molstar_style(**arguments))
    from .controller import manager_for
    from .registry import STYLES, canonical, reference
    from .themes import ANNOTATION_THEMES

    manager = manager_for(_self)
    manager.maintenance()
    style = canonical(style)
    try:
        if style == "list":
            vocabulary = {
                "styles": STYLES,
                "color_themes": (
                    "auto",
                    "keep",
                    *reference()["color_themes"],
                    *ANNOTATION_THEMES,
                ),
                "size_themes": reference()["size_themes"],
            }
            for key, values in vocabulary.items():
                print(f"{key}: " + ", ".join(values))
            for entry in manager.entries.values():
                print(
                    f"{entry.name}: {entry.options['style']}; {max(len(ds) for ds in entry.drawings.values())} states, {entry.seconds:.2f} s, {entry.nbytes / 1024**2:.1f} MiB"
                )
            return vocabulary
        if style == "help":
            print(molstar_style.__doc__)
            return molstar_style.__doc__
        if style in ("reset", "refresh"):
            return getattr(manager, style)(str(name))
        if style in ("png", "ray"):
            from .export import image

            return image(manager, str(filename), width, height, style == "ray")
        entry = manager.apply(
            style,
            str(selection),
            representation,
            str(color),
            str(quality),
            str(name),
            None if transparency in ("keep", "", None) else float(transparency),
            int(state),
            data,
            params,
            float(cache_mb),
            float(gpu_cache_mb),
        )
        if not int(quiet):
            print(
                f"molstar_style: {name} ({style}), prepared in {entry.seconds:.2f} s; {entry.nbytes / 1024**2:.1f} MiB geometry"
            )
        return entry
    except (ValueError, RuntimeError, KeyError, IndexError) as exc:
        from pymol import CmdException

        raise CmdException(str(exc)) from exc


def __init_plugin__(app=None):
    from pymol import cmd
    from pymol.shortcut import Shortcut

    from .controller import manager_for
    from .registry import OPERATIONS, STYLES, reference

    cmd.extend("molstar_style", molstar_style)
    cmd.auto_arg[0]["molstar_style"] = [
        Shortcut((*STYLES, *OPERATIONS)),
        "style or operation",
        "",
    ]
    cmd.auto_arg[1]["molstar_style"] = [cmd.selection_sc, "selection", ""]
    cmd.auto_arg[2]["molstar_style"] = [Shortcut(STYLES), "representation", ""]
    cmd.auto_arg[3]["molstar_style"] = [
        Shortcut(("auto", "keep", *reference()["color_themes"])),
        "color theme",
        "",
    ]
    manager_for(cmd)
