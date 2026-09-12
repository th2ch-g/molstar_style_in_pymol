"""Shared native ray geometry and transactional PNG export."""

from contextlib import contextmanager
from pathlib import Path
from tempfile import NamedTemporaryFile
from uuid import uuid4

import numpy as np

from .mesh import unit


def bake(piece, profile, rotation=None):
    mesh = piece.mesh
    if profile.ignore_light or piece.unlit:
        return mesh.colors
    from .lighting import shade

    rotation = np.eye(3) if rotation is None else rotation
    normals = mesh.normals @ rotation[:3, :3].T
    view = np.broadcast_to([0, 0, 1], normals.shape)
    if rotation.shape == (4, 4):
        view = -(mesh.vertices @ rotation[:3, :3].T + rotation[:3, 3])
    base = shade(mesh.colors, normals, view, profile)
    if profile.params.get("xray"):
        base = mesh.colors * (0.15 + 0.85 * (1 - np.abs(normals[:, 2]))[:, None] ** 1.5)
    return np.clip(base, 0, 1)


def cgo_mesh(piece, profile, rotation=None, ray_only=False):
    from pymol.cgo import ALPHA, BEGIN, COLOR, END, NORMAL, TRIANGLE, TRIANGLES, VERTEX

    mesh = piece.mesh
    # PyMOL reverses BEGIN/TRIANGLES internally; the direct TRIANGLE opcode does not.
    faces = mesh.faces[:, ::-1] if ray_only else mesh.faces
    color = bake(piece, profile, rotation)
    if ray_only:
        values = np.empty((len(faces), 30), np.float32)
        values[:, 0] = ALPHA
        values[:, 1] = mesh.opacity * mesh.alphas[faces].mean(axis=1)
        values[:, 2] = TRIANGLE
        values[:, 3:12] = mesh.vertices[faces].reshape(-1, 9)
        values[:, 12:21] = mesh.normals[faces].reshape(-1, 9)
        values[:, 21:30] = color[faces].reshape(-1, 9)
        return values.ravel().tolist()
    # ALPHA inside BEGIN is not a portable per-vertex ray attribute in PyMOL.
    # Quantized triangle groups retain actual transparency in both CGO paths.
    alpha = np.round(mesh.opacity * mesh.alphas[faces].mean(axis=1) * 255) / 255
    result = []
    for value in np.unique(alpha):
        if value <= 0:
            continue
        ids = faces[alpha == value].ravel()
        values = np.empty((len(ids), 12), np.float32)
        values[:, 0] = NORMAL
        values[:, 1:4] = mesh.normals[ids]
        values[:, 4] = COLOR
        values[:, 5:8] = color[ids]
        values[:, 8] = VERTEX
        values[:, 9:12] = mesh.vertices[ids]
        result.extend(
            [ALPHA, float(value), BEGIN, TRIANGLES, *values.ravel().tolist(), END]
        )
    return result


def ray_proxy(drawing, rotation=None):
    result = []
    for piece in drawing.pieces:
        if opaque(piece):
            result.extend(cgo_mesh(piece, drawing.profile, rotation, ray_only=True))
    from .scene import Piece
    from .volume import density_projection

    for volume in drawing.volumes:
        result.extend(
            cgo_mesh(
                Piece(density_projection(volume), unlit=True),
                drawing.profile,
                ray_only=True,
            )
        )
    return result


def opaque(piece):
    return (
        piece.mesh.opacity >= 0.999999
        and np.min(piece.mesh.alphas, initial=1) >= 0.999999
    )


def ray_cgo(drawing, cmd, include_volumes=True):
    from pymol.cgo import ALPHA, CYLINDER

    from .scene import Piece
    from .volume import density_slices

    matrix = view_matrix(cmd)
    result = []
    for piece in drawing.pieces:
        result.extend(cgo_mesh(piece, drawing.profile, matrix))
        if drawing.profile.edges != "none":
            edges = visible_edges(
                piece.edges, matrix, cmd.get_setting_int("orthoscopic"), False
            )
            values = np.empty((len(edges), 14))
            values[:, 0] = CYLINDER
            values[:, 1:7] = edges[:, :6]
            values[:, 7] = drawing.profile.edge_width / 2
            values[:, 8:11] = drawing.edge_color
            values[:, 11:14] = drawing.edge_color
            result.extend([ALPHA, 1, *values.ravel().tolist()])
    for volume in drawing.volumes if include_volumes else []:
        result.extend(
            cgo_mesh(
                Piece(density_slices(volume, "high", matrix[2, :3]), unlit=True),
                drawing.profile,
                ray_only=True,
            )
        )
    return result


@contextmanager
def settings(cmd, values):
    saved = {key: cmd.get_setting_tuple(key)[1] for key in values}
    try:
        for key, value in values.items():
            cmd.set(key, value)
        yield
    finally:
        for key, value in saved.items():
            cmd.set(key, value if len(value) > 1 else value[0])


def view_matrix(cmd):
    view = np.asarray(cmd.get_view())
    matrix = np.eye(4)
    matrix[:3, :3] = view[:9].reshape(3, 3).T
    matrix[:3, 3] = view[9:12] - matrix[:3, :3] @ view[12:15]
    return matrix


def visible_edges(edges, modelview, orthoscopic, creases):
    if not len(edges):
        return edges
    rotation = modelview[:3, :3]
    mid = 0.5 * (edges[:, :3] + edges[:, 3:6]) @ rotation.T + modelview[:3, 3]
    direction = np.tile([0.0, 0.0, 1.0], (len(mid), 1)) if orthoscopic else unit(-mid)
    a, b = edges[:, 6:9] @ rotation.T, edges[:, 9:12] @ rotation.T
    fa, fb = np.sum(a * direction, axis=1), np.sum(b * direction, axis=1)
    keep = fa * fb <= 0
    if creases:
        keep |= (np.sum(a * b, axis=1) < 0.5) & (np.maximum(fa, fb) > 0)
    return edges[keep]


def has_unmanaged_geometry(cmd, generated):
    """Preserve native lighting when another visible object shares the scene."""
    for name in cmd.get_names("objects", enabled_only=1):
        if name in generated:
            continue
        kind = cmd.get_type(name)
        if kind == "object:molecule":
            rows = []
            cmd.iterate("%" + name, "out.append(reps)", space={"out": rows})
            if any(rows):
                return True
        elif kind not in ("object:group", "object:map", "object:selection"):
            return True
    return False


def image(manager, filename, width, height, ray):
    if not filename:
        raise ValueError("filename is required for PNG and ray export")
    width, height = int(width), int(height)
    if width < 0 or height < 0 or width > 32768 or height > 32768:
        raise ValueError("Image dimensions must be between 0 and 32768 pixels")
    path = Path(filename).expanduser()
    if path.suffix.lower() != ".png":
        path = Path(str(path) + ".png")
    if not path.parent.is_dir():
        raise ValueError("The output directory does not exist")
    cmd = manager.cmd
    playing = cmd.get_movie_playing()
    sculpting = cmd.get_setting_int("sculpting")
    frame = cmd.get_frame()
    show_selection = manager.pool.show_selection
    manager.pool.show_selection = False
    temporary, disabled = [], []
    try:
        active = list(manager.active_drawings())
        panel_list = [p for entry in manager.entries.values() for p in entry.panels]
        has_geometry = any(d.pieces or d.volumes for d in manager.active_drawings())
        if playing:
            cmd.mstop()
        if panel_list and not has_geometry:
            from io import BytesIO

            from .panels import raster

            stream = BytesIO()
            raster(panel_list[0], width or height or 800).save(stream, format="PNG")
            data = stream.getvalue()
        elif ray:
            generated = {name for e in manager.entries.values() for name in e.generated}
            pixel_volumes = any(d.volumes for d in active) and (
                not has_unmanaged_geometry(cmd, generated)
                and all(opaque(p) for d in active for p in d.pieces)
            )
            for drawing in active:
                name = "_molstar_ray_" + uuid4().hex
                payload = ray_cgo(drawing, cmd, include_volumes=not pixel_volumes)
                if payload:
                    temporary.append(name)
                    cmd.load_cgo(payload, name, zoom=0)
            for entry in manager.entries.values():
                for name in entry.generated:
                    if name in cmd.get_names("objects", enabled_only=1):
                        disabled.append(name)
                        cmd.disable(name)
            # Edges are explicit geometry, avoiding an extra global outline pass.
            values = {"ray_trace_mode": 0, "two_sided_lighting": 1}
            if not has_unmanaged_geometry(cmd, set(temporary) | set(disabled)):
                # The material is already evaluated; avoid lighting it a second time.
                values.update(
                    ambient=1,
                    direct=0,
                    reflect=0,
                    specular=0,
                    ray_shadows=0,
                    depth_cue=0,
                    ray_trace_fog=0,
                    ambient_occlusion_mode=0,
                    ray_transparency_specular=0,
                )
            if any(d.volumes for d in active):
                values.update(
                    triangle_max_passes=512,
                    ray_max_passes=1024,
                    ray_transparency_contrast=1,
                    transparency_mode=1,
                    backface_cull=0,
                )
            if any(d.profile.background for d in active):
                values["ray_opaque_background"] = 0
                values["opaque_background"] = 0
            with settings(cmd, values):
                data = cmd.png(None, width, height, ray=1, quiet=1)
            if pixel_volumes:
                from .volume_ray import composite

                data = composite(data, active, cmd)
            from .ray_effects import process

            data = process(data, active, cmd)
        else:
            if manager.widget is None:
                raise ValueError(
                    "GPU PNG export requires the PyMOL Qt GUI; use molstar_style ray in headless mode"
                )
            cmd.draw(width, height, quiet=1)
            data = cmd.png(None, prior=1, quiet=1)
            errors = [
                d.error
                for e in manager.entries.values()
                for ds in e.drawings.values()
                for d in ds
                if d.error
            ]
            if errors:
                raise RuntimeError(errors[0])
        if not isinstance(data, bytes) or not data.startswith(b"\x89PNG\r\n\x1a\n"):
            raise RuntimeError("PyMOL did not return a PNG image")
        with NamedTemporaryFile(dir=path.parent, suffix=".png", delete=False) as stream:
            scratch = Path(stream.name)
            try:
                stream.write(data)
                stream.flush()
                scratch.replace(path)
            finally:
                scratch.unlink(missing_ok=True)
    finally:
        manager.pool.show_selection = show_selection
        for name in temporary:
            cmd.delete(name)
        for name in disabled:
            cmd.enable(name)
        if cmd.get_frame() != frame:
            cmd.frame(frame)
        if sculpting:
            cmd.set("sculpting", sculpting)
        if playing:
            cmd.mplay()
        cmd.rebuild()
        cmd.refresh()
    return str(path)
