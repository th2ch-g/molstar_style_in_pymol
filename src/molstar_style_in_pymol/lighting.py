"""Mol* GGX/Schlick materials and camera-space lights (MIT; see NOTICE)."""

import numpy as np

from .mesh import unit


def settings(profile):
    values = profile.params.get("lighting", {})
    lights = values.get(
        "light", [{"inclination": 150, "azimuth": 320, "intensity": 0.6}]
    )
    if len(lights) > 8:
        raise ValueError("lighting.light supports up to eight lights")
    directions, colors = [], []
    for light in lights:
        inclination = np.deg2rad(float(light.get("inclination", 150)))
        azimuth = np.deg2rad(float(light.get("azimuth", 320)))
        # Mol* negates both the view-space normal and the view direction.
        directions.append(
            -np.array(
                [
                    np.cos(azimuth) * np.sin(inclination),
                    np.sin(azimuth) * np.sin(inclination),
                    np.cos(inclination),
                ]
            )
        )
        colors.append(
            rgb(light.get("color", 0xFFFFFF)) * float(light.get("intensity", 0.6))
        )
    ambient = rgb(values.get("ambientColor", 0xFFFFFF)) * float(
        values.get("ambientIntensity", 0.4)
    )
    exposure = float(values.get("exposure", 1))
    if not np.isfinite(
        np.r_[
            np.asarray(directions).ravel(),
            np.asarray(colors).ravel(),
            ambient,
            exposure,
        ]
    ).all() or np.any(np.r_[np.asarray(colors).ravel(), ambient, exposure] < 0):
        raise ValueError(
            "Lighting angles must be finite; colors and intensities must be finite and nonnegative"
        )
    return (
        np.array(directions).reshape(-1, 3),
        np.array(colors).reshape(-1, 3),
        ambient,
        exposure,
    )


def rgb(value):
    if isinstance(value, (int, np.integer)):
        return np.array([(value >> 16) & 255, (value >> 8) & 255, value & 255]) / 255
    if isinstance(value, str):
        return rgb(int(value.lstrip("#"), 16))
    result = np.asarray(value, float)
    if result.shape != (3,):
        raise ValueError("Lighting color must be an RGB triple or a packed RGB integer")
    return result


def shade(colors, normals, view, profile):
    """Evaluate the same BRDF for native vertex baking and numerical checks."""
    colors, normals, view = np.asarray(colors), unit(normals), unit(view)
    normals = np.where(
        (np.sum(normals * view, axis=-1) < 0)[..., None], -normals, normals
    )
    metal = float(profile.material["metalness"])
    rough = max(0.0525, float(profile.material["roughness"]))
    directions, light_colors, ambient, exposure = settings(profile)
    diffuse = colors * (1 - metal)
    f0 = (1 - metal) * 0.04 + metal * colors
    nv = np.clip(np.sum(normals * view, axis=-1), 0, 1)
    result = diffuse * ambient
    for light, color in zip(directions, light_colors):
        half = unit(light + view)
        nl = np.clip(normals @ light, 0, 1)
        nh = np.clip(np.sum(normals * half, axis=-1), 0, 1)
        vh = np.clip(np.sum(view * half, axis=-1), 0, 1)
        fresnel = np.exp2((-5.55473 * vh - 6.98316) * vh)
        f = f0 * (1 - fresnel[..., None]) + fresnel[..., None]
        a2 = rough**4
        gv = nl * np.sqrt(a2 + (1 - a2) * nv**2)
        gl = nv * np.sqrt(a2 + (1 - a2) * nl**2)
        visibility = 0.5 / np.maximum(gv + gl, 1e-6)
        distribution = a2 / (np.pi * (nh**2 * (a2 - 1) + 1) ** 2)
        specular = f * (visibility * distribution)[..., None]
        if profile.params.get("cel"):
            intensity = nl / np.pi * (1 - metal) + (nl[..., None] * specular) @ [
                0.2126,
                0.7152,
                0.0722,
            ]
            steps = float(profile.params.get("celSteps", 5))
            result += (
                colors * color * np.pi * (np.ceil(intensity * steps) / steps)[..., None]
            )
        else:
            result += color * nl[..., None] * (diffuse + np.pi * specular)
    if metal and not profile.params.get("cel"):
        r = rough * np.array([-1, -0.0275, -0.572, 0.022]) + [1, 0.0425, 1.04, -0.04]
        a004 = np.minimum(r[0] ** 2, np.exp2(-9.28 * nv)) * r[0] + r[1]
        fab = a004[..., None] * [-1.04, 1.04] + r[2:]
        single = f0 * fab[..., :1] + fab[..., 1:]
        ems = 1 - np.sum(fab, axis=-1, keepdims=True)
        average = f0 + (1 - f0) / 21
        multi = single * average / (1 - ems * average) * ems
        irradiance = ambient * metal / np.pi
        result += (
            ambient * metal * single
            + multi * irradiance
            + diffuse * (1 - single - multi) * irradiance
        )
    return np.clip(result, 0.01, 0.99) * exposure
