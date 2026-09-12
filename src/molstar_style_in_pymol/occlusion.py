"""CPU counterpart of the Mol* view-space SSAO used by headless ray export."""

import numpy as np
from scipy.ndimage import map_coordinates

from .mesh import unit
from .postprocessing import ssao_samples


def calculate(distance, cmd, options):
    options = options.get("params", options) if isinstance(options, dict) else {}
    h, w = distance.shape
    valid = np.isfinite(distance)
    if not valid.any():
        return np.ones_like(distance)
    view = cmd.get_view()
    near, far = view[15:17]
    ortho = bool(cmd.get_setting_int("orthoscopic"))
    tangent = np.tan(np.deg2rad(cmd.get_setting_float("field_of_view")) / 2)
    # Use OpenGL's bottom-up convention for the same noise and hemisphere frames.
    distance = distance[::-1].astype(np.float32)
    valid = valid[::-1]
    y, x = np.mgrid[:h, :w].astype(np.float32)
    u, v = (x + 0.5) / w, (y + 0.5) / h
    z = np.where(valid, distance, far)
    half = np.full_like(z, abs(view[11]) * tangent) if ortho else z * tangent
    positions = np.stack(((2 * u - 1) * half * w / h, (2 * v - 1) * half, -z), axis=2)
    depth = (
        (z - near) / (far - near) if ortho else (far - near * far / z) / (far - near)
    )

    def offset(array, dx, dy):
        return array[
            np.clip(np.arange(h) + dy, 0, h - 1)[:, None],
            np.clip(np.arange(w) + dx, 0, w - 1)[None, :],
        ]

    left, right = (
        positions - offset(positions, -1, 0),
        offset(positions, 1, 0) - positions,
    )
    d, up = positions - offset(positions, 0, -1), offset(positions, 0, 1) - positions
    he = abs(2 * offset(depth, -1, 0) - offset(depth, -2, 0) - depth) < abs(
        2 * offset(depth, 1, 0) - offset(depth, 2, 0) - depth
    )
    ve = abs(2 * offset(depth, 0, -1) - offset(depth, 0, -2) - depth) < abs(
        2 * offset(depth, 0, 1) - offset(depth, 0, 2) - depth
    )
    normal = unit(
        np.cross(np.where(he[..., None], left, right), np.where(ve[..., None], d, up))
    ).astype(np.float32)

    def noise(a, b):
        value = np.sin(np.mod(a * 12.9898 + b * 78.233, np.pi)) * 43758.5453
        return value - np.floor(value)

    random = unit(
        np.stack(
            (
                noise(u, v) * 2 - 1,
                noise(u + np.pi, v + 2.71828) * 2 - 1,
                np.zeros_like(u),
            ),
            axis=2,
        )
    ).astype(np.float32)
    t = unit(random - normal * np.sum(random * normal, axis=2)[..., None]).astype(
        np.float32
    )
    b = np.cross(normal, t)
    total, count = np.zeros_like(z), np.zeros_like(z)
    radius = 2 ** float(options.get("radius", 5))
    for sample in ssao_samples():
        point = positions + radius * (
            sample[0] * t + sample[1] * b + sample[2] * normal
        )
        divisor = half if ortho else -point[..., 2] * tangent
        su = point[..., 0] / divisor * h / w * 0.5 + 0.5
        sv = point[..., 1] / divisor * 0.5 + 0.5
        inside = (su >= 0) & (su <= 1) & (sv >= 0) & (sv <= 1)
        coords = [sv * h - 0.5, su * w - 0.5]
        sample_z = map_coordinates(z, coords, order=0, mode="nearest")
        sample_valid = (
            map_coordinates(valid.astype(np.uint8), coords, order=0, mode="nearest") > 0
        )
        ratio = np.clip(radius / np.maximum(abs(z - sample_z), 1e-6), 0, 1)
        smooth = ratio**3 * (ratio * (ratio * 6 - 15) + 10)
        total += inside * sample_valid * (-sample_z >= point[..., 2] + 0.025) * smooth
        count += inside
    result = np.clip(
        1 - float(options.get("bias", 0.8)) * total / np.maximum(count, 1), 0.01, 1
    )
    result[~valid] = 1
    kernel_size = min(25, max(1, int(options.get("blurKernelSize", 15)) | 1))
    pixel_size = half * 2 / h
    depth_bias = float(options.get("blurDepthBias", 0.5))
    for dx, dy in ((1, 0), (0, 1)):
        total, weight = np.zeros_like(z), np.zeros_like(z)
        for i in range(-kernel_size // 2 + 1, kernel_size // 2 + 1):
            sample_z = offset(z, dx * i, dy * i)
            accept = (
                offset(valid, dx * i, dy * i) & valid & (abs(z - sample_z) < depth_bias)
            )
            if abs(i) > 1:
                accept &= abs(i) * pixel_size <= 0.8
            kernel = np.exp(-i * i / (2 * (kernel_size / 3) ** 2))
            weight += kernel * accept
            total += kernel * accept * offset(result, dx * i, dy * i)
        result = np.divide(total, weight, out=np.ones_like(total), where=weight > 0)
    return result[::-1]
