"""Local image, cube-map, and gradient backgrounds shared by export and GPU."""

import numpy as np
from PIL import Image, ImageEnhance, ImageFilter

from .data import local_path
from .themes import rgb


def prepare(value, base="."):
    if value is None:
        return None
    variant = value.get("variant", value)
    name = variant.get("name", "off")
    p = variant.get("params", variant)
    if name == "off":
        return None
    if name not in ("image", "skybox", "horizontalGradient", "radialGradient"):
        raise ValueError("Unknown background variant")

    def read(path):
        image = Image.open(local_path(path, base)).convert("RGB")
        image.thumbnail((2048, 2048))
        if p.get("blur"):
            image = image.filter(ImageFilter.GaussianBlur(float(p["blur"]) * 32))
        image = ImageEnhance.Color(image).enhance(
            max(0, 1 + float(p.get("saturation", 0)))
        )
        image = ImageEnhance.Brightness(image).enhance(
            max(0, 1 + float(p.get("lightness", 0)))
        )
        return np.asarray(image).copy()

    data = {"name": name, "params": p}
    if name == "image":
        src = p.get("source", p.get("file"))
        if isinstance(src, dict):
            src = src.get("params")
        data["image"] = read(src)
    elif name == "skybox":
        faces = p["faces"]
        faces = faces.get("params", faces)
        data["faces"] = {
            k: read(faces[k]) for k in ("px", "nx", "py", "ny", "pz", "nz")
        }
    return data


def raster(data, width, height, rotation=None):
    p = data["params"]
    name = data["name"]
    if name == "image":
        from PIL import ImageOps

        return ImageOps.fit(Image.fromarray(data["image"]), (width, height))
    x, y = np.meshgrid(np.linspace(-1, 1, width), np.linspace(-1, 1, height))
    if name == "skybox":
        directions = np.stack([x * width / height, -y, -np.ones_like(x)], axis=-1)
        if rotation is not None:
            directions = directions @ rotation
        if "rotation" in p:
            from scipy.spatial.transform import Rotation

            angles = p["rotation"]
            directions = (
                directions
                @ Rotation.from_euler(
                    "xyz", [angles.get(k, 0) for k in "xyz"], degrees=True
                ).as_matrix()
            )
        axis = np.argmax(abs(directions), axis=-1)
        out = np.zeros((height, width, 3), np.uint8)
        for k, (positive, negative) in enumerate(
            (("px", "nx"), ("py", "ny"), ("pz", "nz"))
        ):
            for sign, name in ((1, positive), (-1, negative)):
                mask = (axis == k) & (directions[:, :, k] * sign >= 0)
                v = directions[mask]
                denom = np.maximum(abs(v[:, k]), 1e-8)
                if k == 0:
                    u = -sign * v[:, 2] / denom
                    w = -v[:, 1] / denom
                elif k == 1:
                    u = v[:, 0] / denom
                    w = sign * v[:, 2] / denom
                else:
                    u = sign * v[:, 0] / denom
                    w = -v[:, 1] / denom
                image = data["faces"][name]
                ix = np.clip(
                    ((u + 1) / 2 * (image.shape[1] - 1)).astype(int),
                    0,
                    image.shape[1] - 1,
                )
                iy = np.clip(
                    ((w + 1) / 2 * (image.shape[0] - 1)).astype(int),
                    0,
                    image.shape[0] - 1,
                )
                out[mask] = image[iy, ix]
        return Image.fromarray(out)
    ratio = float(p.get("ratio", 0.5))
    if name == "horizontalGradient":
        a, b = rgb(p.get("topColor", 0xDDDDDD)), rgb(p.get("bottomColor", 0xEEEEEE))
        t = (y + 1) / 2
    else:
        a, b = rgb(p.get("centerColor", 0xDDDDDD)), rgb(p.get("edgeColor", 0xEEEEEE))
        t = np.sqrt(x * x + y * y) / np.sqrt(2)
    t = np.clip(t + (0.5 - ratio), 0, 1)
    return Image.fromarray(
        np.uint8((a * (1 - t[:, :, None]) + b * t[:, :, None]) * 255)
    )


def draw(data, pool, projection, modelview, viewport):
    from OpenGL import GL as gl

    # A prepared screen image is cached; cube-map view changes refresh outside draw.
    pixels = data.get("_pixels")
    if pixels is None:
        return
    texture = int(gl.glGenTextures(1))
    try:
        gl.glActiveTexture(gl.GL_TEXTURE0)
        gl.glBindTexture(gl.GL_TEXTURE_2D, texture)
        gl.glTexParameteri(gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MIN_FILTER, gl.GL_LINEAR)
        gl.glTexParameteri(gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MAG_FILTER, gl.GL_LINEAR)
        gl.glTexImage2D(
            gl.GL_TEXTURE_2D,
            0,
            gl.GL_RGB,
            pixels.shape[1],
            pixels.shape[0],
            0,
            gl.GL_RGB,
            gl.GL_UNSIGNED_BYTE,
            np.ascontiguousarray(pixels[::-1]),
        )
        program = pool.program("background")
        gl.glUseProgram(program)
        gl.glUniform1i(gl.glGetUniformLocation(program, "image"), 0)
        gl.glEnable(gl.GL_DEPTH_TEST)
        gl.glDepthFunc(gl.GL_LEQUAL)
        gl.glDepthMask(False)
        gl.glDisable(gl.GL_BLEND)
        gl.glBegin(gl.GL_QUADS)
        for x, y in ((0, 0), (1, 0), (1, 1), (0, 1)):
            gl.glTexCoord2f(x, y)
            gl.glVertex3f(x * 2 - 1, y * 2 - 1, 1)
        gl.glEnd()
    finally:
        gl.glDeleteTextures([texture])
