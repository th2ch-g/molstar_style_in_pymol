"""GLSL scalar-field ray marching with local transfer functions."""

import numpy as np


def draw_volumes(pool, drawing):
    from OpenGL import GL as gl

    modelview, projection, _viewport = drawing.matrices
    matrix = projection @ modelview
    program = pool.program("volume")
    gl.glUseProgram(program)
    for key, value in (("inverseMvp", np.linalg.inv(matrix)), ("mvp", matrix)):
        gl.glUniformMatrix4fv(
            gl.glGetUniformLocation(program, key),
            1,
            True,
            np.asarray(value, np.float32),
        )
    textures = []
    try:
        for volume in drawing.volumes:
            if volume.grid.nbytes > pool.budget:
                raise RuntimeError("Volume texture exceeds gpu_cache_mb")
            cache_key = id(volume)
            if cache_key not in pool.volume_buffers:
                size = volume.pixels.nbytes + volume.lookup.nbytes
                pool.reserve(size)
                pool.volume_buffers[cache_key] = (
                    [int(gl.glGenTextures(1)) for _ in range(2)],
                    size,
                )
                pool.bytes += size
                upload = True
            else:
                upload = False
            pool.volume_buffers.move_to_end(cache_key)
            handles = list(pool.volume_buffers[cache_key][0]) + [
                int(gl.glGenTextures(1))
            ]
            textures.append(handles[2])
            grid = volume.grid
            gl.glActiveTexture(gl.GL_TEXTURE0)
            gl.glBindTexture(gl.GL_TEXTURE_3D, handles[0])
            for key in (gl.GL_TEXTURE_MIN_FILTER, gl.GL_TEXTURE_MAG_FILTER):
                gl.glTexParameteri(gl.GL_TEXTURE_3D, key, gl.GL_LINEAR)
            for key in (
                gl.GL_TEXTURE_WRAP_S,
                gl.GL_TEXTURE_WRAP_T,
                gl.GL_TEXTURE_WRAP_R,
            ):
                gl.glTexParameteri(gl.GL_TEXTURE_3D, key, gl.GL_CLAMP_TO_EDGE)
            if upload:
                gl.glTexImage3D(
                    gl.GL_TEXTURE_3D,
                    0,
                    0x8818,  # GL_LUMINANCE32F_ARB, available with ARB_texture_float.
                    *grid.values.shape,
                    0,
                    gl.GL_LUMINANCE,
                    gl.GL_FLOAT,
                    volume.pixels,
                )
            gl.glUniform1i(gl.glGetUniformLocation(program, "field"), 0)
            gl.glActiveTexture(gl.GL_TEXTURE1)
            gl.glBindTexture(gl.GL_TEXTURE_1D, handles[1])
            lo, hi = volume.transfer[0, 0], volume.transfer[-1, 0]
            gl.glTexParameteri(gl.GL_TEXTURE_1D, gl.GL_TEXTURE_MIN_FILTER, gl.GL_LINEAR)
            gl.glTexParameteri(gl.GL_TEXTURE_1D, gl.GL_TEXTURE_MAG_FILTER, gl.GL_LINEAR)
            gl.glTexParameteri(
                gl.GL_TEXTURE_1D, gl.GL_TEXTURE_WRAP_S, gl.GL_CLAMP_TO_EDGE
            )
            if upload:
                gl.glTexImage1D(
                    gl.GL_TEXTURE_1D,
                    0,
                    gl.GL_RGBA,
                    1024,
                    0,
                    gl.GL_RGBA,
                    gl.GL_FLOAT,
                    volume.lookup,
                )
            gl.glUniform1i(gl.glGetUniformLocation(program, "transfer"), 1)
            gl.glActiveTexture(gl.GL_TEXTURE2)
            gl.glBindTexture(gl.GL_TEXTURE_2D, handles[2])
            current = gl.glGetIntegerv(gl.GL_VIEWPORT)
            gl.glCopyTexImage2D(
                gl.GL_TEXTURE_2D, 0, gl.GL_DEPTH_COMPONENT24, *map(int, current), 0
            )
            gl.glTexParameteri(
                gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MIN_FILTER, gl.GL_NEAREST
            )
            gl.glTexParameteri(
                gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MAG_FILTER, gl.GL_NEAREST
            )
            gl.glUniform1i(gl.glGetUniformLocation(program, "sceneDepth"), 2)
            gl.glUniform3f(
                gl.glGetUniformLocation(program, "dimensions"), *grid.values.shape
            )
            gl.glUniform2f(
                gl.glGetUniformLocation(program, "domain"),
                float(lo),
                float(hi if hi > lo else lo + 1),
            )
            gl.glUniform1f(
                gl.glGetUniformLocation(program, "stepSize"), float(volume.step)
            )
            gl.glUniform1f(
                gl.glGetUniformLocation(program, "opacity"), float(volume.opacity)
            )
            gl.glUniformMatrix4fv(
                gl.glGetUniformLocation(program, "worldToGrid"),
                1,
                True,
                np.asarray(np.linalg.inv(grid.transform), np.float32),
            )
            gl.glActiveTexture(gl.GL_TEXTURE0)
            gl.glEnable(gl.GL_BLEND)
            gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
            gl.glBegin(gl.GL_QUADS)
            for x, y in ((0, 0), (1, 0), (1, 1), (0, 1)):
                gl.glTexCoord2f(x, y)
                gl.glVertex3f(x * 2 - 1, y * 2 - 1, 0)
            gl.glEnd()
    finally:
        if textures:
            gl.glDeleteTextures(textures)
