"""Shared geometric primitives adapted from cuemol_style_in_pymol (MIT)."""

from functools import lru_cache

import numpy as np

from .mesh import Mesh, merge, unit


def indexed(vertices, faces, color, owner=0, opacity=1.0):
    vertices = np.asarray(vertices, float).reshape(-1, 3)
    faces = np.asarray(faces, int).reshape(-1, 3)
    if len(faces) and (faces.min() < 0 or faces.max() >= len(vertices)):
        raise ValueError("Mesh face index is outside the vertex array")
    normals = np.zeros_like(vertices)
    if len(faces):
        p = vertices[faces]
        face_normals = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        for i in range(3):
            np.add.at(normals, faces[:, i], face_normals)
    color = np.broadcast_to(color, (len(vertices), 3))
    owners = np.broadcast_to(owner, (len(vertices),))
    return Mesh(vertices, unit(normals), color, faces, owners, opacity)


def convex(vertices, color, owner=0):
    from scipy.spatial import ConvexHull

    points = np.asarray(vertices, float)
    hull = ConvexHull(points)
    faces = hull.simplices.copy()
    center = points.mean(axis=0)
    for face in faces:
        a, b, c = points[face]
        if np.dot(np.cross(b - a, c - a), a - center) < 0:
            face[1], face[2] = face[2], face[1]
    # Duplicate faces to preserve sharp polyhedral normals.
    return indexed(
        points[faces].reshape(-1, 3), np.arange(faces.size).reshape(-1, 3), color, owner
    )


def cylinder(a, b, radius, color, owner=0, detail=12, radius_b=None):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if np.linalg.norm(b - a) < 1e-9:
        return sphere(a, radius, color, owner, detail)
    axis = unit(b - a)
    x = unit(np.cross(axis, np.eye(3)[np.argmin(abs(axis))]))
    return sweep(
        np.array([a, b]),
        [radius, radius if radius_b is None else radius_b],
        [radius, radius if radius_b is None else radius_b],
        np.array([x, x]),
        np.tile(color, (2, 1)),
        [owner, owner],
        "ellipse",
        detail,
    )


def dashed(a, b, radius, color, owner=0, detail=8, count=7):
    a, b = np.asarray(a), np.asarray(b)
    return merge(
        [
            cylinder(
                a + (b - a) * i / count,
                a + (b - a) * (i + 0.5) / count,
                radius,
                color,
                owner,
                detail,
            )
            for i in range(count)
        ]
    )


def arrow(a, b, radius, color, owner=0, detail=12):
    a, b = np.asarray(a), np.asarray(b)
    length = np.linalg.norm(b - a)
    stem = b - unit(b - a) * min(length * 0.3, 4 * radius)
    return merge(
        [
            cylinder(a, stem, radius, color, owner, detail),
            cylinder(stem, b, radius * 2.2, color, owner, detail, 0.0001),
        ]
    )


def text_mesh(text, position, color, size=0.6, owner=0, axes=None, font=None):
    """Build portable vector text from PyMOL's bundled font."""
    from pymol.vfont import plain

    position = np.asarray(position, float).copy()
    x, y = np.eye(3)[:2] * size if axes is None else np.asarray(axes)[:2] * size
    if any(char not in plain for char in str(text)):
        return unicode_text(str(text), position, color, owner, x, y, font)
    pieces = []
    for char in str(text):
        if char not in plain:
            raise ValueError(
                f"The bundled vector font cannot draw {char!r}; use an ASCII label"
            )
        width, strokes = plain[char]
        last = None
        for move, px, py in np.array(strokes).reshape(-1, 3):
            point = position + x * px + y * py
            if move and last is not None:
                pieces.append(cylinder(last, point, size * 0.035, color, owner, 6))
            last = point
        position += width * x
    return merge(pieces)


def unicode_text(text, position, color, owner, x, y, font=None):
    from PIL import Image, ImageDraw, ImageFont

    if font:
        from .data import local_path

        face = ImageFont.truetype(str(local_path(font)), 48)
        bounds = face.getbbox(text)
        image = Image.new(
            "L", (max(1, bounds[2] - bounds[0] + 4), max(1, bounds[3] - bounds[1] + 4))
        )
        ImageDraw.Draw(image).text(
            (2 - bounds[0], 2 - bounds[1]), text, font=face, fill=255
        )
    else:
        from pymol.Qt import QtCore, QtGui

        if QtGui.QGuiApplication.instance() is None:
            raise ValueError(
                "Unicode labels in headless mode require params.font with a local TTF/OTF file"
            )
        face = QtGui.QFont()
        face.setPixelSize(48)
        metrics = QtGui.QFontMetrics(face)
        width = max(1, metrics.horizontalAdvance(text) + 4)
        height = metrics.height() + 4
        bitmap = QtGui.QImage(width, height, QtGui.QImage.Format_ARGB32)
        bitmap.fill(0)
        painter = QtGui.QPainter(bitmap)
        painter.setFont(face)
        painter.setPen(QtGui.QColor("white"))
        painter.drawText(2, metrics.ascent() + 2, text)
        painter.end()
        buffer = QtCore.QBuffer()
        buffer.open(QtCore.QIODevice.WriteOnly)
        bitmap.save(buffer, "PNG")
        from io import BytesIO

        image = Image.open(BytesIO(bytes(buffer.data()))).getchannel("A")
    mask = np.asarray(image) > 127
    vertices = []
    faces = []
    for row in range(mask.shape[0]):
        edges = np.diff(np.r_[False, mask[row], False].astype(int))
        starts = np.flatnonzero(edges == 1)
        ends = np.flatnonzero(edges == -1)
        for start, end in zip(starts, ends):
            base = len(vertices)
            vertices.extend(
                position + x * col / 48 + y * (mask.shape[0] - r) / 48
                for col, r in (
                    (start, row),
                    (end, row),
                    (end, row + 1),
                    (start, row + 1),
                )
            )
            faces.extend(((base, base + 1, base + 2), (base, base + 2, base + 3)))
    return indexed(vertices, faces, color, owner)


@lru_cache(maxsize=32)
def section(kind, detail):
    if kind == "rectangle":
        # Duplicate corners to keep face normals sharp.
        points = np.array(
            [[-1, -1], [1, -1], [1, -1], [1, 1], [1, 1], [-1, 1], [-1, 1], [-1, -1]],
            float,
        )
        normals = np.array(
            [[0, -1], [0, -1], [1, 0], [1, 0], [0, 1], [0, 1], [-1, 0], [-1, 0]], float
        )
        return points, normals
    if kind == "fancy":
        # Two circular rails joined by flat faces, using Fancy1's sharp=0.3.
        theta = 0.3 * np.pi
        angle = np.linspace(
            theta - np.pi / 2, 3 * np.pi / 2 - theta, max(8, detail // 2) + 1
        )
        rail = np.c_[(1.1 + 0.2 * np.sin(angle)) / 1.3, np.cos(angle)]
        normal = np.c_[1.3 * np.sin(angle), 0.2 * np.cos(angle)]
        points = np.vstack((rail, rail[-1], -rail[0], -rail, -rail[-1], rail[0]))[::-1]
        normals = np.vstack((normal, [0, -1], [0, -1], -normal, [0, 1], [0, 1]))[::-1]
        return points, unit(normals)
    t = np.linspace(0, 2 * np.pi, detail, endpoint=False)
    points = np.c_[np.cos(t), np.sin(t)]
    return points, points.copy()


def frames(path, hints, tangent=None):
    tangent = unit(np.gradient(path, axis=0) if tangent is None else tangent)
    side = hints - tangent * np.sum(hints * tangent, axis=1, keepdims=True)
    for i in range(len(side)):
        if np.linalg.norm(side[i]) < 1e-7:
            previous = side[i - 1] if i else np.eye(3)[np.argmin(np.abs(tangent[i]))]
            side[i] = previous - tangent[i] * np.dot(previous, tangent[i])
        if i and np.dot(side[i], side[i - 1]) < 0:
            side[i] *= -1
    side = unit(side)
    return side, unit(np.cross(tangent, side))


def sweep(
    path,
    widths,
    thickness,
    hints,
    colors,
    owners,
    kind,
    detail,
    back=False,
    front=None,
    side_color=False,
    frame=None,
):
    if len(path) < 2:
        return Mesh([], [], [], [], [])
    side, up = frames(path, hints) if frame is None else frame
    shape, sn = section(kind, detail)
    k = len(shape)
    width = np.broadcast_to(widths, (len(path),))
    thick = np.broadcast_to(thickness, (len(path),))
    v = (
        path[:, None]
        + side[:, None] * width[:, None, None] * shape[None, :, 0, None]
        + up[:, None] * thick[:, None, None] * shape[None, :, 1, None]
    )
    n = unit(
        side[:, None] * sn[None, :, 0, None] / np.maximum(width[:, None, None], 1e-6)
        + up[:, None] * sn[None, :, 1, None] / np.maximum(thick[:, None, None], 1e-6)
    )
    col = np.repeat(np.asarray(colors)[:, None, :], k, axis=1)
    if back or side_color:
        # Frame signs can depend on preceding strands and loops. Color helix
        # undersides by the physical inward direction without twisting the mesh.
        polarity = (
            np.where(np.sum(up * front, axis=1) < 0, -1, 1)
            if front is not None
            else np.ones(len(path))
        )
        mask = polarity[:, None] * shape[None, :, 1] < -0.1
        if front is not None:
            mask &= np.linalg.norm(front, axis=1)[:, None] > 1e-7
        if side_color:
            mask = np.broadcast_to(np.abs(sn[:, 0]) > 0.5, col.shape[:2])
        original = col[mask]
        value = original.max(axis=-1, keepdims=True)
        saturation = (value - original.min(axis=-1, keepdims=True)) / np.maximum(
            value, 1e-8
        )
        scale = np.maximum(saturation - 0.4, 0) / np.maximum(saturation, 1e-8)
        col[mask] = value - (value - original) * scale
    j = np.flatnonzero(
        np.linalg.norm(shape - np.roll(shape, -1, axis=0), axis=1) > 1e-8
    )
    a = (np.arange(len(path) - 1)[:, None] * k + j).ravel()
    b = (np.arange(len(path) - 1)[:, None] * k + (j + 1) % k).ravel()
    f = np.concatenate((np.c_[a, b, a + k], np.c_[b, b + k, a + k])).tolist()
    vertices, normals, vertex_colors = (
        v.reshape(-1, 3),
        n.reshape(-1, 3),
        col.reshape(-1, 3),
    )
    ids = np.repeat(owners, k)
    # Independent cap vertices prevent shading the cross-section as a side wall.
    for idx, sign in ((0, -1), (-1, 1)):
        rim = v[idx]
        start = len(vertices)
        normal = unit(np.cross(side[idx], up[idx]))
        vertices = np.vstack((vertices, path[idx], rim))
        normals = np.vstack((normals, np.tile(sign * normal, (k + 1, 1))))
        vertex_colors = np.vstack((vertex_colors, np.tile(colors[idx], (k + 1, 1))))
        ids = np.r_[ids, np.full(k + 1, owners[idx])]
        for j in range(k):
            tri = (start, start + 1 + j, start + 1 + (j + 1) % k)
            f.append(tri if sign > 0 else tri[::-1])
    faces = np.asarray(f)
    # Correct winding against analytic normals, including capped section seams.
    fn = np.cross(
        vertices[faces[:, 1]] - vertices[faces[:, 0]],
        vertices[faces[:, 2]] - vertices[faces[:, 0]],
    )
    reverse = np.sum(fn * normals[faces].mean(axis=1), axis=1) < 0
    faces[reverse] = faces[reverse, ::-1]
    return Mesh(vertices, normals, vertex_colors, faces, ids)


@lru_cache(maxsize=8)
def sphere_template(detail):
    # Rings exclude the poles, whose triangles are added separately.
    phi = np.linspace(0, np.pi, detail // 2 + 1)[1:-1]
    theta = np.linspace(0, 2 * np.pi, detail, endpoint=False)
    v = np.c_[
        np.outer(np.sin(phi), np.cos(theta)).ravel(),
        np.outer(np.sin(phi), np.sin(theta)).ravel(),
        np.repeat(np.cos(phi), detail),
    ]
    f = []
    for i in range(len(phi) - 1):
        for j in range(detail):
            a = i * detail + j
            b = i * detail + (j + 1) % detail
            f.extend(((a, b, a + detail), (b, b + detail, a + detail)))
    top, bottom = len(v), len(v) + 1
    for j in range(detail):
        f.extend(
            (
                (top, (j + 1) % detail, j),
                (
                    bottom,
                    (len(phi) - 1) * detail + j,
                    (len(phi) - 1) * detail + (j + 1) % detail,
                ),
            )
        )
    v = np.vstack((v, [0, 0, 1], [0, 0, -1]))
    f = np.asarray(f)
    reverse = (
        np.sum(
            np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
            * v[f].mean(axis=1),
            axis=1,
        )
        < 0
    )
    f[reverse] = f[reverse, ::-1]
    return v, f


def sphere(center, radius, color, owner, detail):
    v, f = sphere_template(detail)
    return Mesh(
        v * radius + center, v, np.tile(color, (len(v), 1)), f, np.full(len(v), owner)
    )


def bond(a, b, radius, colors, owners, detail):
    if np.linalg.norm(b - a) < 1e-8:
        return Mesh([], [], [], [], [])
    t = unit(b - a)
    h = np.eye(3)[np.argmin(np.abs(t))]
    middle = (a + b) / 2
    # Independent half cylinders keep the element-color boundary sharp.
    return merge(
        [
            sweep(
                np.array([start, end]),
                radius,
                radius,
                np.tile(h, (2, 1)),
                [color, color],
                [owner, owner],
                "ellipse",
                detail,
            )
            for start, end, color, owner in (
                (a, middle, colors[0], owners[0]),
                (middle, b, colors[1], owners[1]),
            )
        ]
    )
