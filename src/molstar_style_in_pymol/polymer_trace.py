"""Mol* residue-local polymer curves and cross sections (MIT; see NOTICE)."""

import numpy as np

from .mesh import Mesh, unit
from .primitives import sphere


def orthogonal(tangent, direction):
    result = direction - tangent * np.dot(tangent, direction)
    if np.linalg.norm(result) < 1e-10:
        axis = np.array([1.0, 0, 0])
        if abs(tangent[0]) > 0.999:
            axis = np.array([0.0, 1, 0])
        result = axis - tangent * np.dot(tangent, axis)
    return unit(result)


def match(vector, reference):
    return -vector if np.dot(vector, reference) < 0 else vector


def spline(points, t, tension):
    a, b, c, d = points
    v0, v1 = (c - a) * tension, (d - b) * tension
    return (
        (2 * b - 2 * c + v0 + v1) * t**3
        + (-3 * b + 3 * c - 2 * v0 - v1) * t**2
        + v0 * t
        + b
    )


def interpolate(controls, tension=0.5, shift=0.5, linear=8):
    """Port curve-segment.ts, including boundary tensions and normal smoothing."""
    points = np.array([controls[f"p{i}"] for i in range(5)])
    begin = 0.5 if controls["secStrucFirst"] else tension
    end = 0.5 if controls["secStrucLast"] else tension
    curve, tangents = [], []
    for t in np.linspace(0, 1, linear + 1):
        if t < 1 - shift:
            p, s = points[:4], t + shift
            te = begin + (tension - begin) * t
            tt = begin
        else:
            p, s = points[1:], t - (1 - shift)
            te = tension + (end - tension) * t
            tt = te
        curve.append(spline(p, s, te))
        tangents.append(unit(spline(p, s + 0.01, tt) - spline(p, s - 0.01, tt)))
    curve = np.asarray(curve, np.float32)
    tangents = np.asarray(tangents, np.float32)
    first = orthogonal(tangents[0], np.asarray(controls["d12"]))
    last = match(orthogonal(tangents[-1], np.asarray(controls["d23"])), first)
    previous = first
    normals = np.zeros_like(curve)
    binormals = np.zeros_like(curve)
    for i, tangent in enumerate(tangents):
        fraction = i / linear
        j = fraction * fraction * (3 - 2 * fraction) * linear
        t = 0 if i == 0 else 1 / (linear + 1 - j)
        dot = np.clip(np.dot(previous, last), -1, 1)
        theta = np.arccos(dot) * t
        direction = previous * np.cos(theta) + unit(last - previous * dot) * np.sin(
            theta
        )
        normal = orthogonal(tangent, direction)
        normals[i] = normal
        binormals[i] = unit(np.cross(tangent, normal))
        previous = normal
    for i in range(1, linear):
        normals[i] = (normals[i - 1] + normals[i] + normals[i + 1]) / 3
        binormals[i] = unit(np.cross(tangents[i], normals[i]))
    return curve, normals, binormals


def helix_centers(points):
    """GROMACS helixorient centers, as used by Mol* tubular helices."""
    if len(points) < 4:
        return points.copy()
    centers = np.zeros_like(points)
    for i in range(len(points) - 3):
        a, b, c, d = points[i : i + 4]
        x, y = 2 * b - a - c, 2 * c - b - d
        lx, ly = np.linalg.norm(x), np.linalg.norm(y)
        if min(lx, ly) < 1e-10:
            centers[i + 1 : i + 3] = [b, c]
            continue
        radius = np.sqrt(lx * ly) / max(2, 2 * (1 - np.dot(x, y) / (lx * ly)))
        centers[i + 1] = b - x * radius / lx
        centers[i + 2] = c - y * radius / ly
    for index, adjacent, next_index in ((0, 1, 2), (-1, -2, -3)):
        axis = unit(centers[adjacent] - centers[next_index])
        centers[index] = centers[adjacent] + axis * np.dot(
            points[index] - centers[adjacent], axis
        )
    return centers


def controls_for(points, directions, secondary, index, tubular=False):
    """Port the atomic trace iterator for one continuous, noncyclic chain."""
    n = len(points)
    indices = np.clip(np.arange(index - 3, index + 4), 0, n - 1)
    ss = np.array(secondary)[indices]
    p = points[indices].copy()
    helix = secondary[index] == "H"
    if tubular:
        # Keep center-line segments from pulling neighboring loops into a helix.
        for side in (-1, 1):
            for distance in (1, 2, 3):
                k = 3 + side * distance
                if (ss[k] == "H") != helix:
                    if helix:
                        for j in range(distance, 4):
                            p[3 + side * j] = p[3 + side * (distance - 1)]
                    else:
                        delta = 2 * (p[k] - p[k - side])
                        for j in range(distance, 4):
                            p[3 + side * j] = p[3 + side * (j - 1)] + delta
                    break
    dirs = directions[np.clip(np.arange(index - 1, index + 3), 0, n - 1)].copy()
    if np.linalg.norm(directions[index]) < 1e-10:
        for k in range(4):
            dirs[k] = unit(np.cross(p[k + 2] - p[k + 1], p[k + 3] - p[k + 1]))
    # Extend control points past termini by 1.5 Å, before sheet smoothing.
    for side, available in ((-1, index), (1, n - index - 1)):
        for distance in (1, 2, 3):
            k = 3 + side * distance
            boundary = distance > available or (tubular and helix and ss[k] != ss[3])
            if boundary:
                delta = unit(p[k - side] - p[k - 2 * side]) * 1.5
                for j in range(distance, 4):
                    p[3 + side * j] = p[3 + side * (j - 1)] + delta
                break
    controls = {}
    for k in range(1, 6):
        controls[f"p{k - 1}"] = (
            (p[k - 1] + 2 * p[k] + p[k + 1]) / 4
            if ss[k] == "S" or (tubular and ss[k] == "H")
            else p[k]
        )
    for k, key in ((1, "d12"), (2, "d23")):
        controls[key] = (
            match(dirs[k - 1], dirs[k]) + 2 * dirs[k] + match(dirs[k + 1], dirs[k])
        ) / 4
    controls.update(
        secStrucFirst=index > 0 and secondary[index - 1] != secondary[index],
        secStrucLast=index == n - 1
        or secondary[min(index + 1, n - 1)] != secondary[index],
        first=index == 0,
        last=index == n - 1,
        initial=index == 0,
        final=index == n - 1,
    )
    return controls


def section_mesh(
    curve,
    normals,
    binormals,
    widths,
    heights,
    color,
    owner,
    profile,
    radial,
    caps,
    arrow=0,
    round_cap=False,
):
    """Use Mol*'s section dimensions, opposite diagonals and sharp sheet faces."""
    count = len(curve)
    cap_factor = None
    if round_cap and any(caps):
        index = np.arange(count, dtype=float)
        span = count - 1
        if all(caps):
            cap_start = index <= span / 2
            fraction = np.abs(index - span / 2) / (span / 2)
        else:
            cap_start = np.full(count, caps[0])
            fraction = (span - index if caps[0] else index) / span
        cap_factor = np.sqrt(np.maximum(0, 1 - fraction**2))
        widths, heights = widths * cap_factor, heights * cap_factor
    square = profile == "square" or radial == 4
    if square or radial == 2:
        section = np.array(
            [[1, 1], [1, -1], [1, -1], [-1, -1], [-1, -1], [-1, 1], [-1, 1], [1, 1]],
            float,
        )
        sn = np.array(
            [[1, 0], [1, 0], [0, -1], [0, -1], [-1, 0], [-1, 0], [0, 1], [0, 1]], float
        )
        if radial == 2:
            section[:, 1] = 0
            sn = np.tile([0, 1], (8, 1))
        h = arrow * np.linspace(1, 0, count) if arrow else heights
        vertices = (
            curve[:, None]
            + normals[:, None] * h[:, None, None] * section[None, :, 0, None]
            + binormals[:, None] * widths[:, None, None] * section[None, :, 1, None]
        )
        ns = (
            normals[:, None] * sn[None, :, 0, None]
            + binormals[:, None] * sn[None, :, 1, None]
        )
        if arrow:
            slope = arrow / max(np.linalg.norm(curve[-1] - curve[0]), 1e-10)
            ns[:, [0, 1, 4, 5]] += unit(np.cross(normals, binormals))[:, None] * slope
        edges = [0, 2, 4, 6]
    else:
        angle = (np.arange(radial) * 2 + (profile == "rounded")) * np.pi / radial
        cosine, sine = np.cos(angle), np.sin(angle)
        vertices = (
            curve[:, None]
            + normals[:, None] * heights[:, None, None] * cosine[None, :, None]
            + binormals[:, None] * widths[:, None, None] * sine[None, :, None]
        )
        ns = (
            normals[:, None] * widths[:, None, None] * cosine[None, :, None]
            + binormals[:, None] * heights[:, None, None] * sine[None, :, None]
        )
        if profile == "rounded":
            offset = np.maximum(heights - widths, 0)
            vertices = (
                curve[:, None]
                + normals[:, None]
                * (widths[:, None] * cosine + offset[:, None] * np.sign(cosine))[
                    ..., None
                ]
                + binormals[:, None] * widths[:, None, None] * sine[None, :, None]
            )
            ns = (
                normals[:, None] * cosine[None, :, None]
                + binormals[:, None] * sine[None, :, None]
            )
            q1 = int(np.floor(radial / 4 + 0.5))
            for j, sign in ((q1 - 1, 1), (q1, 1), (3 * q1 - 1, -1), (3 * q1, -1)):
                if j < radial:
                    ns[:, j] = sign * binormals
        edges = np.arange(radial)
    ns = unit(ns)
    if cap_factor is not None:
        cap_normal = (
            unit(np.cross(normals, binormals)) * np.where(cap_start, -1, 1)[:, None]
        )
        dot = np.clip(np.sum(cap_normal[:, None] * ns, axis=2), -1, 1)
        theta = np.arccos(dot) * cap_factor[:, None]
        relative = unit(ns - cap_normal[:, None] * dot[..., None])
        ns = (
            cap_normal[:, None] * np.cos(theta)[..., None]
            + relative * np.sin(theta)[..., None]
        )
    k = vertices.shape[1]
    faces = []
    for i in range(count - 1):
        for j in edges:
            a, b = i * k + j, i * k + (j + 1) % k
            if j < k / 2:
                faces.extend(((a, b, b + k), (a, b + k, a + k)))
            else:
                faces.extend(((a, b, a + k), (b, b + k, a + k)))
    vertices = vertices.reshape(-1, 3)
    ns = ns.reshape(-1, 3)
    for index, enabled, sign in ((0, caps[0], -1), (count - 1, caps[1], 1)):
        if not enabled or radial == 2 or (square and arrow and sign == 1):
            continue
        base = len(vertices)
        if square:
            vertices = np.vstack(
                (vertices, vertices[index * k + np.array([0, 1, 3, 5])])
            )
            ns = np.vstack(
                (
                    ns,
                    np.tile(
                        sign * unit(np.cross(normals[index], binormals[index])), (4, 1)
                    ),
                )
            )
            faces.extend(((base, base + 1, base + 2), (base, base + 2, base + 3)))
            continue
        vertices = np.vstack(
            (vertices, curve[index], vertices[index * k : (index + 1) * k])
        )
        ns = np.vstack(
            (
                ns,
                np.tile(
                    sign * unit(np.cross(normals[index], binormals[index])), (k + 1, 1)
                ),
            )
        )
        faces.extend((base, base + 1 + j, base + 1 + (j + 1) % k) for j in edges)
    if square and arrow and not caps[0]:
        # Close both shoulders of the enlarged arrow without capping the ribbon.
        for sign in (-1, 1):
            base = len(vertices)
            shoulder = np.array(
                [
                    curve[0] + sign * normals[0] * h + side * binormals[0] * widths[0]
                    for h, side in (
                        (heights[0], -1),
                        (arrow, -1),
                        (arrow, 1),
                        (heights[0], 1),
                    )
                ]
            )
            vertices = np.vstack((vertices, shoulder))
            ns = np.vstack(
                (ns, np.tile(-unit(np.cross(normals[0], binormals[0])), (4, 1)))
            )
            faces.extend(((base, base + 1, base + 2), (base, base + 2, base + 3)))
    faces = np.array(faces)
    # PyMOL and Mol* use opposite normal signs in their fragment-lighting code.
    face_normals = np.cross(
        vertices[faces[:, 1]] - vertices[faces[:, 0]],
        vertices[faces[:, 2]] - vertices[faces[:, 0]],
    )
    reverse = np.sum(face_normals * ns[faces].mean(axis=1), axis=1) < 0
    faces[reverse] = faces[reverse, ::-1]
    return Mesh(
        vertices,
        ns,
        np.tile(color, (len(vertices), 1)),
        faces,
        np.full(len(vertices), owner),
    )


def trace_meshes(
    state, ids, group_lookup, colors, radii, params, linear, radial, putty=False
):
    points = state.coords[ids].copy()
    atoms = [state.atoms[i] for i in ids]
    secondary = ["" if putty else a.ss for a in atoms]
    directions = np.zeros_like(points)
    for k, i in enumerate(ids):
        names = {state.atoms[j].name: j for j in group_lookup[i]}
        pair = ("C", "O")
        if atoms[k].kind == "nucleic":
            pair = ("C3'", "C1'") if atoms[k].resn.startswith("D") else ("C4'", "C3'")
        if all(n in names for n in pair):
            directions[k] = state.coords[names[pair[1]]] - state.coords[names[pair[0]]]
    tubular = bool(params.get("tubularHelices", False)) and not putty
    if tubular:
        centers = helix_centers(points)
        for k, ss in enumerate(secondary):
            if ss == "H":
                points[k] = centers[k]
                directions[k] = [1, 0, 0]
    size = float(params.get("sizeFactor", 0.2))
    aspect = float(params.get("aspectRatio", 5))
    arrow_factor = float(params.get("arrowFactor", 1.5))
    for key in ("helixProfile", "nucleicProfile"):
        if params.get(key, "elliptical") not in ("elliptical", "rounded", "square"):
            raise ValueError(f"{key} must be elliptical, rounded or square")
    for index, atom in enumerate(atoms):
        control = controls_for(points, directions, secondary, index, tubular)
        nucleic = atom.kind == "nucleic"
        shift = 0.3 if nucleic else 0.5
        curve, normals, binormals = interpolate(
            control,
            0.9 if atom.ss == "H" and not (tubular or putty) else 0.5,
            shift,
            linear,
        )
        neighbors = np.clip([index - 1, index, index + 1], 0, len(ids) - 1)
        widths = size * radii[np.array(ids)[neighbors]]
        if putty and "bfactorScale" in params:
            widths = size * np.array(
                [
                    max(
                        0.25,
                        np.sqrt(
                            max(atoms[k].bfactor, 0) / float(params["bfactorScale"])
                        ),
                    )
                    for k in neighbors
                ]
            )
        if len(ids) == 1:
            yield sphere(points[0], widths[1] * 2, colors[ids[0]], ids[0], radial)
            continue
        helix = atom.ss == "H" and not putty
        sheet = atom.ss == "S" and not putty
        heights = (
            widths * aspect
            if (helix or sheet or nucleic) and not putty
            else widths.copy()
        )
        if helix and tubular:
            heights = widths = widths * aspect * 1.5
        t = np.linspace(0, 1, linear + 1)
        samples = np.where(t < 1 - shift, t + shift, t - (1 - shift) + 1)
        width = np.interp(samples, [0, 1, 2], widths)
        height = np.interp(samples, [0, 1, 2], heights)
        segment_count = linear
        if index == 0:
            segment_count = max(int(np.floor(linear * shift + 0.5)), 1)
            offset = linear - segment_count
            curve, normals, binormals = (
                curve[offset:].copy(),
                normals[offset:].copy(),
                binormals[offset:].copy(),
            )
            curve[0] = control["p2"] + unit(control["p2"] - curve[1]) * widths[1] * 2
        elif index == len(ids) - 1:
            segment_count = max(int(np.floor(linear * (1 - shift) + 0.5)), 1)
            curve = curve[: segment_count + 1].copy()
            normals, binormals = (
                normals[: segment_count + 1],
                binormals[: segment_count + 1],
            )
            curve[-1] = control["p2"] + unit(control["p2"] - curve[-2]) * widths[1] * 2
        # Upstream keeps the first size samples after trimming terminal curves.
        width, height = width[: segment_count + 1], height[: segment_count + 1]
        profile = (
            "square"
            if sheet
            else params.get(
                "nucleicProfile" if nucleic else "helixProfile",
                "square" if nucleic else "elliptical",
            )
        )
        if np.allclose(width, height):
            profile = "elliptical"
        if nucleic:
            normals, binormals = -binormals, normals
        caps = (
            control["first"] or control["secStrucFirst"],
            control["last"] or control["secStrucLast"],
        )
        arrow = heights[1] * arrow_factor if sheet and control["secStrucLast"] else 0
        yield section_mesh(
            curve,
            normals,
            binormals,
            width,
            height,
            colors[ids[index]],
            ids[index],
            profile,
            radial,
            caps,
            arrow,
            round_cap=helix and tubular and bool(params.get("roundCap", False)),
        )
