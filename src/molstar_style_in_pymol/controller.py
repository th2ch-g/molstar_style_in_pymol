"""Transactional named views, shared source visibility, and session recovery."""

import re
from dataclasses import dataclass, field
from time import perf_counter
from uuid import uuid4

import numpy as np

from . import export, geometry, source
from .data import parameters, pymol_grid, read_data
from .formats import hydrate
from .gpu import Drawing, Pool
from .registry import PRESETS, QUALITIES, STRUCTURE, resolve
from .scene import PickTarget


@dataclass
class Entry:
    name: str
    options: dict
    saved: dict
    drawings: dict
    sources: tuple
    generated: list = field(default_factory=list)
    panels: list = field(default_factory=list)
    windows: list = field(default_factory=list)
    seconds: float = 0
    native_bytes: int = 0
    source_view: object = None

    @property
    def nbytes(self):
        return sum(
            p.mesh.nbytes + p.edges.nbytes
            for ds in self.drawings.values()
            for d in ds
            for p in d.pieces
        ) + sum(
            v.grid.nbytes
            for ds in self.drawings.values()
            for d in ds
            for v in d.volumes
        )


def combine(states):
    from copy import deepcopy

    from chempy.models import Indexed

    model = Indexed()
    atoms = []
    coords = []
    reps = []
    bonds = []
    for state in states:
        offset = len(atoms)
        atoms.extend(state.atoms)
        coords.extend(state.coords)
        reps.extend(state.representations)
        local = deepcopy(state.model)
        model.atom.extend(local.atom)
        for bond in local.bond:
            bond.index = [int(i) + offset for i in bond.index]
            model.bond.append(bond)
            bonds.append(tuple(bond.index))
    return source.State(
        tuple(atoms),
        np.asarray(coords, np.float32).reshape(-1, 3),
        tuple(bonds),
        model,
        np.asarray(reps),
        {},
    )


class Manager:
    def __init__(self, cmd):
        self.cmd = cmd
        self.entries = {}
        self.saved = {}
        self.pool = Pool()
        self.widget = self.events = None
        self.busy = False
        self.selected_keys = None

    def attach(self):
        if self.events is None:
            from .picking import attach

            self.widget, self.events = attach(self)

    def release_gpu(self):
        if self.pool.context is not None and self.widget is not None:
            self.cmd._call_with_opengl_context(self.pool.clear)

    def load_entry(self, entry, group):
        cmd = self.cmd

        def load(values, name, state):
            entry.native_bytes += len(values) * 4
            if (
                entry.native_bytes + entry.nbytes
                > float(entry.options["cache_mb"]) * 1024**2
            ):
                raise ValueError("Geometry and native CGO together exceed cache_mb")
            cmd.load_cgo(values, name, state=state, zoom=0)

        for k, drawings in enumerate(entry.drawings.values()):
            shape = f"{group}_shape_{k}"
            alpha = f"{group}_alpha_{k}"
            ray = f"{group}_ray_{k}"
            has_alpha = any(not export.opaque(p) for d in drawings for p in d.pieces)
            has_ray = any(
                d.volumes or any(export.opaque(p) for p in d.pieces) for d in drawings
            )
            names = (
                [shape] + ([alpha] if has_alpha else []) + ([ray] if has_ray else [])
            )
            entry.generated.extend(names)
            for index, drawing in enumerate(drawings, 1):
                drawing.name = shape
                drawing.state = index
                cmd.load_callback(drawing, shape, index, 1, 0, 1, 0)
                if has_alpha:
                    values = []
                    for piece in drawing.pieces:
                        if not export.opaque(piece):
                            values.extend(
                                export.cgo_mesh(
                                    piece,
                                    drawing.profile,
                                    export.view_matrix(cmd),
                                    ray_only=True,
                                )
                            )
                    load(values, alpha, index)
                if has_ray:
                    load(export.ray_proxy(drawing, export.view_matrix(cmd)), ray, index)
            for name in names:
                cmd.group(group, name)
            if has_alpha:
                cmd.set("cgo_lighting", 0, alpha)
                cmd.set("cgo_transparency", 0, alpha)
        if not entry.generated:
            cmd.group(group)

    def unload_entry(self, entry, delete_group=True):
        for window in entry.windows:
            if window is not None:
                window.close()
        for name in entry.generated:
            self.cmd.delete(name)
        if delete_group:
            self.cmd.delete(entry.name)

    def apply(
        self,
        style,
        selection="all",
        representation=None,
        color="auto",
        quality="medium",
        name="molstar",
        transparency=None,
        state=0,
        data=None,
        params=None,
        cache_mb=2048,
        gpu_cache_mb=256,
    ):
        if (
            not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name)
            or self.cmd.get_legal_name(name) != name
        ):
            raise ValueError(
                "name must be an unreserved identifier starting with a letter"
            )
        if quality not in QUALITIES:
            raise ValueError(f"quality must be one of {tuple(QUALITIES)}")
        if (
            not np.isfinite(float(cache_mb))
            or float(cache_mb) <= 0
            or not np.isfinite(float(gpu_cache_mb))
            or float(gpu_cache_mb) <= 0
        ):
            raise ValueError("Cache budgets must be finite and positive")
        if transparency is not None and (
            not np.isfinite(float(transparency)) or not 0 <= float(transparency) <= 1
        ):
            raise ValueError("transparency must be between zero and one")
        options = {
            "style": style,
            "selection": selection,
            "representation": representation,
            "color": color,
            "quality": quality,
            "name": name,
            "transparency": transparency,
            "state": state,
            "data": data,
            "params": params,
            "cache_mb": cache_mb,
            "gpu_cache_mb": gpu_cache_mb,
        }
        p = parameters(params)
        profile = resolve(style, representation, p)
        planes = np.asarray(profile.params.get("clipPlanes", []), float)
        if planes.size and (
            planes.ndim != 2
            or planes.shape[1] != 4
            or len(planes) > 6
            or not np.isfinite(planes).all()
        ):
            raise ValueError("clipPlanes accepts at most six [nx,ny,nz,offset] planes")
        payload = read_data(data)
        from .background import prepare, raster

        profile.background = prepare(p.get("background"), payload.get("_base", "."))
        if profile.background:
            profile.background["_pixels"] = np.asarray(
                raster(profile.background, 512, 512)
            )
        if (
            not payload
            and selection in self.cmd.get_names("objects")
            and self.cmd.get_type(selection) == "object:map"
        ):
            payload = {
                "grid": pymol_grid(
                    self.cmd, selection, int(state) if int(state) > 0 else 1
                )
            }
        if "tables" in payload and profile.representation in (
            "direct-volume",
            "dot",
            "isosurface",
            "segment",
            "slice",
        ):
            from .data import read_grid

            payload["grid"] = read_grid(payload["file"])
        budget = int(float(cache_mb) * 1024**2)
        profile.params["_budget"] = budget
        started = perf_counter()
        snapshots = {}
        saved = {}
        molecular = profile.representation in (
            *STRUCTURE,
            *PRESETS,
            "interactions",
            "cross-link-restraint",
            "clashes",
        )
        molecular_objects = (
            self.cmd.get_object_list(f"({selection})")
            if selection and profile.representation not in ("mvs", "kinemage", "g3d")
            else []
        )
        if molecular or molecular_objects:
            try:
                snapshots, saved = source.read(self.cmd, selection, budget)
            except ValueError:
                if molecular:
                    raise
        if name in self.cmd.get_names("objects") and name not in self.entries:
            raise ValueError(f"An unrelated object named {name!r} already exists")
        saved = {key: self.saved.get(key, value) for key, value in saved.items()}
        drawings = []
        panel_list = []
        total = 0
        nstates = max((len(v) for v in snapshots.values()), default=1)
        if int(state) < 0 or int(state) > nstates:
            raise ValueError(f"state must be 0 (all) or 1..{nstates}")
        indices = range(nstates) if int(state) == 0 else [int(state) - 1]
        for index in indices:
            snapshot = (
                combine(
                    [
                        values[min(index, len(values) - 1)]
                        for values in snapshots.values()
                    ]
                )
                if snapshots
                else None
            )
            annotations = hydrate(
                {**payload, "state": index + 1}, snapshot.atoms if snapshot else ()
            )
            if (
                profile.representation == "unitcell"
                and "cell" not in annotations
                and snapshots
            ):
                annotations["cell"] = self.cmd.get_symmetry(next(iter(snapshots)))[:6]
            g = geometry.build(
                snapshot,
                profile.representation,
                color,
                profile.params,
                annotations,
                quality,
            )
            if len(planes):
                from .clipping import clip

                for piece in g.pieces:
                    piece.mesh = clip(piece.mesh, planes)
                g.pieces = [piece for piece in g.pieces if len(piece.mesh.faces)]
                for v in g.volumes:
                    from .data import Grid

                    points = v.grid.world(
                        np.indices(v.grid.values.shape).reshape(3, -1).T
                    )
                    mask = np.all(
                        points @ planes[:, :3].T + planes[:, 3] >= 0, axis=1
                    ).reshape(v.grid.values.shape)
                    floor = float(v.grid.values.min() - max(1, np.ptp(v.grid.values)))
                    v.grid = Grid(
                        np.where(mask, v.grid.values, floor),
                        v.grid.transform.copy(),
                        v.grid.label,
                    )
                    v.transfer = np.r_[[[floor, *v.transfer[0, 1:4], 0]], v.transfer]
                    v.prepare()
            for piece in g.pieces:
                if not piece.atoms:
                    piece.atoms = (PickTarget(label=profile.representation),)
                if transparency is not None:
                    piece.mesh.opacity *= 1 - float(transparency)
                if profile.edges != "none":
                    piece.edges = piece.mesh.edge_data()
            for v in g.volumes:
                if transparency is not None:
                    v.opacity *= 1 - float(transparency)
            total += g.nbytes
            if total > budget:
                raise ValueError(
                    "Prepared geometry exceeds cache_mb; reduce quality, selection, or states"
                )
            d = Drawing(g.pieces, profile, (0, 0, 0), self.pool, volumes=g.volumes)
            d.anchors = snapshot.coords if snapshot else np.empty((0, 3), np.float32)
            d.keys = (
                tuple((a.model, a.index) for a in snapshot.atoms) if snapshot else ()
            )
            if g.volumes:
                corners = np.concatenate(
                    [
                        v.grid.world(
                            np.array(
                                [
                                    [x, y, z]
                                    for x in (0, v.grid.values.shape[0] - 1)
                                    for y in (0, v.grid.values.shape[1] - 1)
                                    for z in (0, v.grid.values.shape[2] - 1)
                                ]
                            )
                        )
                        for v in g.volumes
                    ]
                )
                d.extent = [corners.min(axis=0).tolist(), corners.max(axis=0).tolist()]
            drawings.append(d)
            panel_list.extend(g.panels)
        entry = Entry(
            name,
            options,
            saved,
            {"selection": drawings},
            tuple(snapshots),
            panels=panel_list,
        )
        staging = f"molstar_stage_{uuid4().hex}"
        old = self.entries.get(name)
        self.busy = True
        try:
            self.attach()
            self.load_entry(entry, staging)
            if entry.nbytes > budget:
                raise ValueError("Prepared mesh cache exceeds cache_mb")
            self.cmd.disable(staging)
            # Every failure-prone preparation and load finishes before replacement.
            source.hide_reps(self.cmd, saved)
            if old:
                self.unload_entry(old)
            self.cmd.set_name(staging, name)
            self.cmd.enable(name)
            self.entries[name] = entry
            self.saved.update(saved)
            self.restore_unused()
            self.release_gpu()
            self.pool.budget = int(float(gpu_cache_mb) * 1024**2)
            self.selected_keys = None
            from .panels import show

            for panel in panel_list:
                entry.windows.append(show(panel, self.cmd))
        except Exception:
            if self.entries.get(name) is not entry:
                self.unload_entry(entry, delete_group=False)
                self.cmd.delete(staging)
                source.restore_reps(
                    self.cmd, {k: v for k, v in saved.items() if k not in self.saved}
                )
            raise
        finally:
            self.busy = False
        entry.seconds = perf_counter() - started
        self.cmd.refresh()
        return entry

    def restore_unused(self):
        used = (
            set().union(*(e.saved for e in self.entries.values()))
            if self.entries
            else set()
        )
        unused = {k: v for k, v in self.saved.items() if k not in used}
        source.restore_reps(self.cmd, unused)
        for key in unused:
            self.saved.pop(key, None)

    def reset(self, name):
        names = list(self.entries) if name == "all" else [name]
        for key in names:
            if key not in self.entries:
                if name != "all":
                    raise ValueError(f"No managed view named {key!r}")
                continue
            entry = self.entries.pop(key)
            self.unload_entry(entry)
        self.restore_unused()
        self.release_gpu()
        self.cmd.refresh()

    def refresh(self, name):
        for key in list(self.entries) if name == "all" else [name]:
            if key not in self.entries:
                raise ValueError(f"No managed view named {key!r}")
            self.apply(**self.entries[key].options)

    def active_drawings(self):
        visible = self.cmd.get_vis()
        for entry in self.entries.values():
            if not visible.get(entry.name, [0])[0]:
                continue
            for drawings in entry.drawings.values():
                name = drawings[0].name
                if not visible.get(name, [0])[0]:
                    continue
                if self.cmd.get_setting_int("all_states", name):
                    yield from drawings
                else:
                    index = (
                        self.cmd.get_object_state(name) - 1 if len(drawings) > 1 else 0
                    )
                    if 0 <= index < len(drawings):
                        yield drawings[index]

    def maintenance(self):
        if self.busy:
            return
        names = set(self.cmd.get_names("objects"))
        for key, entry in list(self.entries.items()):
            errors = [d.error for ds in entry.drawings.values() for d in ds if d.error]
            if (
                key not in names
                or any(n not in names for n in entry.generated)
                or any(n not in names for n in entry.sources)
                or errors
            ):
                self.reset(key)
                if errors:
                    print(
                        f"molstar_style: restored source after drawing failure: {errors[0]}"
                    )
                continue
            if len(entry.sources) == 1:
                obj = entry.sources[0]
                visible = self.cmd.get_vis().get(obj, [0])[0]
                current = (
                    self.cmd.get_setting_int("state", obj),
                    self.cmd.get_setting_int("all_states", obj),
                    visible,
                )
                if current != entry.source_view:
                    for target in entry.generated:
                        self.cmd.set("state", current[0], target)
                        self.cmd.set("all_states", current[1], target)
                        (self.cmd.enable if visible else self.cmd.disable)(target)
                    entry.source_view = current
        self.update_selection()
        view = tuple(self.cmd.get_view())
        for entry in self.entries.values():
            for ds in entry.drawings.values():
                for d in ds:
                    background = d.profile.background
                    if (
                        background
                        and background["name"] == "skybox"
                        and background.get("_view") != view
                    ):
                        from .background import raster

                        background["_pixels"] = np.asarray(
                            raster(
                                background, 512, 512, np.array(view[:9]).reshape(3, 3)
                            )
                        )
                        background["_view"] = view
                        self.cmd.refresh()

    def update_selection(self):
        names = [
            n
            for n in self.cmd.get_names("selections", enabled_only=1)
            if not n.startswith("_")
        ]
        keys = (
            frozenset(self.cmd.index(" or ".join("%" + n for n in names)))
            if names
            else frozenset()
        )
        if keys == self.selected_keys:
            return
        self.selected_keys = keys
        for entry in self.entries.values():
            for ds in entry.drawings.values():
                for d in ds:
                    indices = [i for i, key in enumerate(d.keys) if key in keys]
                    d.selected = np.ascontiguousarray(d.anchors[indices])
        self.cmd.refresh()


def manager_for(cmd):
    owner = cmd._pymol
    manager = vars(owner).get("_molstar_style_manager")
    if manager is None:
        manager = owner._molstar_style_manager = Manager(cmd)
    for attr, hook in (
        ("_session_save_tasks", save_session),
        ("_session_restore_tasks", restore_session),
    ):
        tasks = getattr(owner, attr)
        if hook not in tasks:
            tasks.append(hook)
    return manager


def save_session(session, _self):
    manager = manager_for(_self)
    session["molstar_style"] = {
        "version": 1,
        "saved": dict(manager.saved),
        "entries": [
            {"options": dict(e.options), "generated": list(e.generated)}
            for e in manager.entries.values()
        ],
    }
    return 1


def restore_session(session, _self):
    manager = manager_for(_self)
    if manager.events is not None:
        manager.events.close()
    manager.release_gpu()
    manager.entries.clear()
    manager.saved.clear()
    manager.widget = manager.events = None
    data = session.get("molstar_style", {})
    source.restore_reps(_self, data.get("saved", {}))
    for entry in data.get("entries", []):
        for name in entry["generated"]:
            _self.delete(name)
        _self.delete(entry["options"]["name"])
    for entry in data.get("entries", []):
        try:
            manager.apply(**entry["options"])
        except Exception as exc:
            print(
                f"molstar_style: source restored; managed view could not be reloaded: {exc}"
            )
    return 1
