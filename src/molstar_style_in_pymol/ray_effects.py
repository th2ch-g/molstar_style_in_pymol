"""Image-space approximations for native-ray export, using mesh depth samples."""

from io import BytesIO

import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter, minimum_filter, shift

from .postprocessing import effect_value


def depth_samples(drawings, cmd, width, height):
    from .export import view_matrix

    scale = min(1, 256 / max(width, height))
    w = max(2, int(width * scale))
    h = max(2, int(height * scale))
    matrix = view_matrix(cmd)
    view = cmd.get_view()
    fov = np.deg2rad(cmd.get_setting_float("field_of_view"))
    half = abs(view[11]) * np.tan(fov / 2)
    depth = np.full((h, w), np.inf)
    orthoscopic = cmd.get_setting_int("orthoscopic")
    for d in drawings:
        for p in d.pieces:
            points = p.mesh.vertices @ matrix[:3, :3].T + matrix[:3, 3]
            divisor = (
                np.full(len(points), half)
                if orthoscopic
                else np.maximum(1e-6, -points[:, 2]) * np.tan(fov / 2)
            )
            x = (points[:, 0] / divisor * height / width + 1) / 2 * (w - 1)
            y = (1 - points[:, 1] / divisor) / 2 * (h - 1)
            projected = np.c_[x, y, -points[:, 2]]
            for triangle in projected[p.mesh.faces]:
                if (triangle[:, 2] <= 0).any():
                    continue
                lo = np.maximum(0, np.ceil(triangle[:, :2].min(axis=0)).astype(int))
                hi = np.minimum(
                    [w - 1, h - 1], np.floor(triangle[:, :2].max(axis=0)).astype(int)
                )
                if np.any(lo > hi):
                    continue
                a, b, c = triangle
                denominator = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (
                    a[1] - c[1]
                )
                if abs(denominator) < 1e-12:
                    continue
                yy, xx = np.mgrid[lo[1] : hi[1] + 1, lo[0] : hi[0] + 1]
                wa = (
                    (b[1] - c[1]) * (xx - c[0]) + (c[0] - b[0]) * (yy - c[1])
                ) / denominator
                wb = (
                    (c[1] - a[1]) * (xx - c[0]) + (a[0] - c[0]) * (yy - c[1])
                ) / denominator
                wc = 1 - wa - wb
                inside = (wa >= -1e-7) & (wb >= -1e-7) & (wc >= -1e-7)
                if orthoscopic:
                    values = wa * a[2] + wb * b[2] + wc * c[2]
                else:
                    inverse = wa / a[2] + wb / b[2] + wc / c[2]
                    values = np.divide(
                        1, inverse, out=np.full_like(inverse, np.inf), where=inverse > 0
                    )
                region = depth[lo[1] : hi[1] + 1, lo[0] : hi[0] + 1]
                np.minimum(region, np.where(inside, values, np.inf), out=region)
    mask = np.isfinite(depth)
    if not mask.any():
        return np.ones((height, width))
    lo, hi = float(depth[mask].min()), float(depth[mask].max())
    filled = np.ones_like(depth)
    filled[mask] = (depth[mask] - lo) / max(hi - lo, 1e-6) * 0.8 + 0.1
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
    original_rgb = rgb.copy()
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
    if not cmd.get_setting_int("bg_gradient"):
        setting_type, value = cmd.get_setting_tuple("bg_rgb")
        background_color = np.asarray(
            cmd.get_color_tuple(value[0]) if setting_type == 5 else value
        )
        untouched = (alpha > 0.999) & (
            np.max(abs(original_rgb - background_color), axis=2) < 1 / 255
        )
        rgb[untouched] = original_rgb[untouched]
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
