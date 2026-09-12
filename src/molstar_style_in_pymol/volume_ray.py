"""Pixel-wise scalar-field integration for dedicated headless ray output."""

from io import BytesIO

import numpy as np
from PIL import Image
from scipy.ndimage import map_coordinates

from .mesh import unit
from .volume import transfer_values


def render(volume, cmd, depth, width, height):
    """Integrate up to the nearest opaque mesh, using the GPU sampling rule."""
    from .export import view_matrix

    matrix = view_matrix(cmd)
    view = cmd.get_view()
    ortho = bool(cmd.get_setting_int("orthoscopic"))
    tangent = np.tan(np.deg2rad(cmd.get_setting_float("field_of_view")) / 2)
    inverse = np.linalg.inv(volume.grid.transform)
    output = np.zeros((height, width, 4), np.float32)
    # Bound working memory independently of the requested image height.
    for top in range(0, height, 32):
        rows = min(32, height - top)
        yy, xx = np.mgrid[top : top + rows, :width]
        x = ((xx.ravel() + 0.5) / width * 2 - 1) * tangent * width / height
        y = (1 - (yy.ravel() + 0.5) / height * 2) * tangent
        if ortho:
            origin = np.c_[x * abs(view[11]), y * abs(view[11]), np.zeros(x.size)]
            direction = np.tile([0, 0, -1], (x.size, 1))
        else:
            origin = np.zeros((x.size, 3))
            direction = unit(np.c_[x, y, -np.ones(x.size)])
        forward = -direction[:, 2]
        origin = (origin - matrix[:3, 3]) @ matrix[:3, :3]
        direction = direction @ matrix[:3, :3]
        o = origin @ inverse[:3, :3].T + inverse[:3, 3]
        d = direction @ inverse[:3, :3].T
        safe = np.where(d >= 0, 1, -1) * np.maximum(abs(d), 1e-9)
        a, b = -o / safe, (np.asarray(volume.grid.values.shape) - 1 - o) / safe
        start = np.maximum(np.minimum(a, b).max(axis=1), view[15] / forward)
        end = np.minimum(np.maximum(a, b).min(axis=1), view[16] / forward)
        end = np.minimum(end, depth[top : top + rows].ravel() / forward)
        delta = np.maximum(volume.step, (end - start) / 1023)
        total = np.zeros((x.size, 4), float)
        for step in range(1024):
            t = start + step * delta
            active = (start < end) & (t <= end) & (total[:, 3] <= 0.985)
            if not active.any():
                break
            indices = o[active] + t[active, None] * d[active]
            values = map_coordinates(
                volume.grid.values, indices.T, order=1, mode="nearest"
            )
            rgba = transfer_values(volume, values, delta[active])
            if volume.color_grid is not None:
                rgba[:, :3] = np.stack(
                    [
                        map_coordinates(
                            volume.color_grid[..., k],
                            indices.T,
                            order=1,
                            mode="nearest",
                        )
                        for k in range(3)
                    ],
                    axis=1,
                )
            contribution = (1 - total[active, 3]) * rgba[:, 3]
            total[active, :3] += rgba[:, :3] * contribution[:, None]
            total[active, 3] += contribution
        total[total[:, 3] < 0.005] = 0
        output[top : top + rows] = total.reshape(rows, width, 4)
    return output


def composite(data, drawings, cmd):
    """Combine integrated density with native mesh ray output and its background."""
    from .ray_effects import depth_samples

    image = Image.open(BytesIO(data)).convert("RGBA")
    width, height = image.size
    array = np.asarray(image, dtype=float) / 255
    array[..., :3] *= array[..., 3:]
    depth = depth_samples(drawings, cmd, width, height, full=True, raw=True)
    for drawing in drawings:
        for volume in drawing.volumes:
            layer = render(volume, cmd, depth, width, height)
            array = layer + array * (1 - layer[..., 3:])
    array[..., :3] = np.divide(
        array[..., :3],
        array[..., 3:],
        out=np.zeros_like(array[..., :3]),
        where=array[..., 3:] > 0,
    )
    stream = BytesIO()
    Image.fromarray(np.uint8(np.clip(array, 0, 1) * 255), "RGBA").save(
        stream, format="PNG"
    )
    return stream.getvalue()
