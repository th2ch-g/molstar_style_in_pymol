"""Clip the shared mesh, so native ray and GPU use identical cut geometry."""

import numpy as np

from .mesh import Mesh, unit


def clip(mesh, planes):
    vertices = []
    normals = []
    colors = []
    alphas = []
    owners = []
    faces = []
    for face in mesh.faces:
        polygon = [
            np.r_[
                mesh.vertices[i],
                mesh.normals[i],
                mesh.colors[i],
                mesh.alphas[i],
                mesh.owners[i],
            ]
            for i in face
        ]
        for plane in planes:
            output = []
            if not polygon:
                break
            for a, b in zip(polygon, polygon[1:] + polygon[:1]):
                da = np.dot(a[:3], plane[:3]) + plane[3]
                db = np.dot(b[:3], plane[:3]) + plane[3]
                if da >= 0:
                    output.append(a)
                if (da >= 0) != (db >= 0):
                    t = da / (da - db)
                    point = a + (b - a) * t
                    point[-1] = a[-1] if t < 0.5 else b[-1]
                    output.append(point)
            polygon = output
        if len(polygon) < 3:
            continue
        base = len(vertices)
        for row in polygon:
            vertices.append(row[:3])
            normals.append(row[3:6])
            colors.append(row[6:9])
            alphas.append(row[9])
            owners.append(int(row[10]))
        faces.extend([base, base + i, base + i + 1] for i in range(1, len(polygon) - 1))
    return Mesh(
        vertices,
        unit(np.array(normals).reshape(-1, 3)),
        colors,
        faces,
        owners,
        mesh.opacity,
        alphas,
    )
