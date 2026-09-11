"""Image-space approximations for native-ray export, using mesh depth samples."""

from io import BytesIO

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, minimum_filter, shift

from .postprocessing import effect_value


def depth_samples(drawings, cmd, width, height):
    from scipy.spatial import cKDTree

    from .export import view_matrix

    scale = min(1, 256 / max(width, height))
    w = max(2, int(width * scale))
    h = max(2, int(height * scale))
    matrix = view_matrix(cmd)
    view = cmd.get_view()
    fov = np.deg2rad(cmd.get_setting_float("field_of_view"))
    half = abs(view[11]) * np.tan(fov / 2)
    samples = []
    for d in drawings:
        for p in d.pieces:
            points = p.mesh.vertices @ matrix[:3, :3].T + matrix[:3, 3]
            divisor = (
                np.full(len(points), half)
                if cmd.get_setting_int("orthoscopic")
                else np.maximum(1e-6, -points[:, 2]) * np.tan(fov / 2)
            )
            x = (points[:, 0] / divisor * height / width + 1) / 2 * (w - 1)
            y = (1 - points[:, 1] / divisor) / 2 * (h - 1)
            samples.append(np.c_[x, y, -points[:, 2]])
    if not samples:
        return np.ones((height, width))
    points = np.concatenate(samples)
    good = (
        (points[:, 0] >= 0)
        & (points[:, 0] < w)
        & (points[:, 1] >= 0)
        & (points[:, 1] < h)
    )
    points = points[good]
    if not len(points):
        return np.ones((height, width))
    depth = np.full((h, w), np.inf)
    indices = np.rint(points[:, :2]).astype(int)
    np.minimum.at(depth, (indices[:, 1], indices[:, 0]), points[:, 2])
    yy, xx = np.indices(depth.shape)
    mask = np.isfinite(depth)
    distance, index = cKDTree(np.c_[xx[mask], yy[mask]]).query(
        np.c_[xx.ravel(), yy.ravel()]
    )
    filled = depth[mask][index].reshape(depth.shape)
    lo, hi = float(points[:, 2].min()), float(points[:, 2].max())
    filled = (filled - lo) / max(hi - lo, 1e-6) * 0.8 + 0.1
    filled[distance.reshape(depth.shape) > 4] = 1
    return np.asarray(
        Image.fromarray(filled.astype(np.float32)).resize(
            (width, height), Image.Resampling.BILINEAR
        )
    )


def process(data, drawings, cmd):
    effects = {}
    background = None
    for d in drawings:
        effects.update(d.profile.effects)
        background = d.profile.background or background
    image = Image.open(BytesIO(data)).convert("RGBA")
    width, height = image.size
    array = np.asarray(image).astype(float) / 255
    rgb = array[:, :, :3]
    alpha = array[:, :, 3]
    depth = depth_samples(drawings, cmd, width, height)
    strength = effect_value(effects.get("occlusion")) + 0.5 * effect_value(
        effects.get("illumination")
    )
    if strength:
        near = minimum_filter(depth, size=max(3, int(width / 100) | 1))
        delta = depth - near
        occlusion = gaussian_filter(
            ((delta > 0.004) & (delta < 0.2)).astype(float), max(0.5, width / 400)
        )
        rgb *= 1 - np.clip(occlusion * strength * 0.3, 0, 0.7)[:, :, None]
    strength = effect_value(effects.get("shadow"))
    if strength:
        other = shift(depth, [height / 80, -width / 100], mode="constant", cval=1)
        rgb *= 1 - ((other < depth - 0.02) & (other > depth - 0.3))[:, :, None] * min(
            0.5, strength * 0.3
        )
    strength = effect_value(effects.get("bloom"))
    if strength:
        bright = rgb * np.maximum(0, rgb.max(axis=2) - 0.65)[:, :, None]
        rgb += (
            gaussian_filter(bright, [max(1, width / 100), max(1, width / 100), 0])
            * strength
        )
    strength = effect_value(effects.get("dof"))
    if strength:
        blurred = gaussian_filter(rgb, [max(1, width / 200), max(1, width / 200), 0])
        weight = np.clip(
            abs(depth - float(effects.get("focus", 0.5))) * strength * 2, 0, 0.9
        )
        rgb = rgb * (1 - weight[:, :, None]) + blurred * weight[:, :, None]
    image = Image.fromarray(
        np.uint8(np.clip(np.dstack([rgb, alpha]), 0, 1) * 255), "RGBA"
    )
    if background:
        from .background import raster

        base = raster(
            background, width, height, np.array(cmd.get_view()[:9]).reshape(3, 3)
        ).convert("RGBA")
        opacity = float(background["params"].get("opacity", 1))
        if opacity < 1:
            from .themes import rgb as color

            base = Image.blend(
                Image.new(
                    "RGBA",
                    (width, height),
                    tuple(np.uint8(color(cmd.get_setting_tuple("bg_rgb")[1]) * 255))
                    + (255,),
                ),
                base,
                opacity,
            )
        base.alpha_composite(image)
        image = base
    stream = BytesIO()
    image.save(stream, format="PNG")
    return stream.getvalue()
