"""Scoped screen effects with explicit Qt framebuffer restoration."""

from functools import lru_cache


def effect_value(value):
    if isinstance(value, dict):
        if value.get("name") == "off" or value.get("enabled") is False:
            return 0.0
        value = value.get("params", value)
        return float(value.get("strength", value.get("intensity", 1)))
    return float(value or 0)


def api():
    from OpenGL import GL as gl

    if gl.glGenFramebuffers:
        return (
            gl.glGenFramebuffers,
            gl.glBindFramebuffer,
            gl.glFramebufferTexture2D,
            gl.glCheckFramebufferStatus,
            gl.glDeleteFramebuffers,
        )
    from OpenGL.GL.EXT import framebuffer_object as ext

    return (
        ext.glGenFramebuffersEXT,
        ext.glBindFramebufferEXT,
        ext.glFramebufferTexture2DEXT,
        ext.glCheckFramebufferStatusEXT,
        ext.glDeleteFramebuffersEXT,
    )


class Framebuffer:
    def __init__(self, width, height):
        from OpenGL import GL as gl

        self.gen, self.bind, self.attach, self.check, self.delete = api()
        self.previous = int(gl.glGetIntegerv(gl.GL_FRAMEBUFFER_BINDING))
        self.viewport = gl.glGetIntegerv(gl.GL_VIEWPORT).copy()
        self.width, self.height = width, height
        self.handle = int(self.gen(1))
        self.textures = []
        self.auxiliary = []
        self.bind(gl.GL_FRAMEBUFFER, self.handle)
        try:
            for attachment, internal, fmt in (
                (gl.GL_COLOR_ATTACHMENT0, gl.GL_RGBA, gl.GL_RGBA),
                (
                    gl.GL_DEPTH_ATTACHMENT,
                    gl.GL_DEPTH_COMPONENT24,
                    gl.GL_DEPTH_COMPONENT,
                ),
            ):
                texture = int(gl.glGenTextures(1))
                self.textures.append(texture)
                gl.glBindTexture(gl.GL_TEXTURE_2D, texture)
                gl.glTexParameteri(
                    gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MIN_FILTER, gl.GL_NEAREST
                )
                gl.glTexParameteri(
                    gl.GL_TEXTURE_2D, gl.GL_TEXTURE_MAG_FILTER, gl.GL_NEAREST
                )
                gl.glTexParameteri(
                    gl.GL_TEXTURE_2D, gl.GL_TEXTURE_WRAP_S, gl.GL_CLAMP_TO_EDGE
                )
                gl.glTexParameteri(
                    gl.GL_TEXTURE_2D, gl.GL_TEXTURE_WRAP_T, gl.GL_CLAMP_TO_EDGE
                )
                gl.glTexImage2D(
                    gl.GL_TEXTURE_2D,
                    0,
                    internal,
                    width,
                    height,
                    0,
                    fmt,
                    gl.GL_FLOAT,
                    None,
                )
                self.attach(gl.GL_FRAMEBUFFER, attachment, gl.GL_TEXTURE_2D, texture, 0)
            if self.check(gl.GL_FRAMEBUFFER) != gl.GL_FRAMEBUFFER_COMPLETE:
                raise RuntimeError("OpenGL framebuffer is incomplete")
        except Exception:
            self.close()
            raise
        self.restore()

    def begin(self):
        from OpenGL import GL as gl

        self.previous = int(gl.glGetIntegerv(gl.GL_FRAMEBUFFER_BINDING))
        self.viewport = gl.glGetIntegerv(gl.GL_VIEWPORT).copy()
        self.bind(gl.GL_FRAMEBUFFER, self.handle)
        gl.glViewport(0, 0, self.width, self.height)
        gl.glClearColor(0, 0, 0, 0)
        gl.glClearDepth(1)
        gl.glDepthMask(True)
        gl.glClear(gl.GL_COLOR_BUFFER_BIT | gl.GL_DEPTH_BUFFER_BIT)

    def ambient_occlusion(self, pool, effects, projection):
        import numpy as np
        from OpenGL import GL as gl

        if not self.auxiliary:
            self.auxiliary = [Framebuffer(self.width, self.height) for _ in range(2)]
        value = effects.get("occlusion", {})
        p = value.get("params", value) if isinstance(value, dict) else {}
        inverse = np.linalg.inv(projection).astype(np.float32)
        source = self.textures[1]
        for index, name in enumerate(("ssao", "ssao_blur", "ssao_blur")):
            target = self.auxiliary[index % 2]
            program = pool.program(name)
            target.begin()
            gl.glUseProgram(program)
            gl.glDisable(gl.GL_DEPTH_TEST)
            gl.glDisable(gl.GL_BLEND)
            for slot, (key, texture) in enumerate(
                (("depth", self.textures[1]), ("ao", source))
            ):
                gl.glActiveTexture(gl.GL_TEXTURE0 + slot)
                gl.glBindTexture(gl.GL_TEXTURE_2D, texture)
                gl.glUniform1i(gl.glGetUniformLocation(program, key), slot)
            gl.glUniformMatrix4fv(
                gl.glGetUniformLocation(program, "inverseProjection"), 1, True, inverse
            )
            gl.glUniform2f(
                gl.glGetUniformLocation(program, "pixel"),
                1 / self.width,
                1 / self.height,
            )
            if index == 0:
                gl.glUniformMatrix4fv(
                    gl.glGetUniformLocation(program, "projection"),
                    1,
                    True,
                    projection.astype(np.float32),
                )
                gl.glUniform1f(
                    gl.glGetUniformLocation(program, "radius"),
                    2 ** float(p.get("radius", 5)),
                )
                gl.glUniform1f(
                    gl.glGetUniformLocation(program, "bias"), float(p.get("bias", 0.8))
                )
                gl.glUniform3fv(
                    gl.glGetUniformLocation(program, "samples"), 32, ssao_samples()
                )
            else:
                gl.glUniform2f(
                    gl.glGetUniformLocation(program, "direction"),
                    int(index == 1),
                    int(index == 2),
                )
                gl.glUniform1f(
                    gl.glGetUniformLocation(program, "depthBias"),
                    float(p.get("blurDepthBias", 0.5)),
                )
                gl.glUniform1i(
                    gl.glGetUniformLocation(program, "kernelSize"),
                    min(25, max(1, int(p.get("blurKernelSize", 15)) | 1)),
                )
            quad()
            source = target.textures[0]
            target.restore()
        return source

    def composite(self, program, effects, projection, pool):
        from OpenGL import GL as gl

        ao = self.textures[0]
        if effect_value(effects.get("occlusion")) > 0:
            ao = self.ambient_occlusion(pool, effects, projection)
        self.bind(gl.GL_FRAMEBUFFER, self.previous)
        gl.glViewport(*self.viewport)
        gl.glUseProgram(program)
        for unit, (name, texture) in enumerate(zip(("image", "depth"), self.textures)):
            gl.glActiveTexture(gl.GL_TEXTURE0 + unit)
            gl.glBindTexture(gl.GL_TEXTURE_2D, texture)
            gl.glUniform1i(gl.glGetUniformLocation(program, name), unit)
        gl.glActiveTexture(gl.GL_TEXTURE2)
        gl.glBindTexture(gl.GL_TEXTURE_2D, ao)
        gl.glUniform1i(gl.glGetUniformLocation(program, "ambientOcclusion"), 2)
        gl.glUniform2f(
            gl.glGetUniformLocation(program, "pixel"), 1 / self.width, 1 / self.height
        )
        for name in (
            "occlusion",
            "outline",
            "shadow",
            "bloom",
            "dof",
            "antialias",
            "illumination",
        ):
            gl.glUniform1f(
                gl.glGetUniformLocation(program, name), effect_value(effects.get(name))
            )
        gl.glUniform1f(
            gl.glGetUniformLocation(program, "focus"), float(effects.get("focus", 0.5))
        )
        gl.glActiveTexture(gl.GL_TEXTURE0)
        gl.glEnable(gl.GL_BLEND)
        gl.glBlendFuncSeparate(
            gl.GL_SRC_ALPHA,
            gl.GL_ONE_MINUS_SRC_ALPHA,
            gl.GL_ONE,
            gl.GL_ONE_MINUS_SRC_ALPHA,
        )
        gl.glEnable(gl.GL_DEPTH_TEST)
        gl.glDepthMask(True)
        quad()

    def close(self):
        from OpenGL import GL as gl

        for target in self.auxiliary:
            target.close()
        self.auxiliary = []
        if int(gl.glGetIntegerv(gl.GL_FRAMEBUFFER_BINDING)) == self.handle:
            self.bind(gl.GL_FRAMEBUFFER, self.previous)
        if self.textures:
            gl.glDeleteTextures(self.textures)
            self.textures = []
        if self.handle:
            self.delete(1, [self.handle])
            self.handle = 0

    def restore(self):
        from OpenGL import GL as gl

        self.bind(gl.GL_FRAMEBUFFER, self.previous)
        gl.glViewport(*self.viewport)


def quad():
    from OpenGL import GL as gl

    gl.glActiveTexture(gl.GL_TEXTURE0)
    gl.glBegin(gl.GL_QUADS)
    for x, y in ((0, 0), (1, 0), (1, 1), (0, 1)):
        gl.glTexCoord2f(x, y)
        gl.glVertex3f(x * 2 - 1, y * 2 - 1, 0)
    gl.glEnd()


@lru_cache(maxsize=1)
def ssao_samples():
    """Mol*'s seeded PCG blue-noise hemisphere and quadratic distance ramp."""
    import numpy as np

    state = 26699

    def random():
        nonlocal state
        old = state
        state = (state * 1664525 + 1013904223) & 0xFFFFFFFF
        x = ((old >> 18) ^ old) >> 5
        rotation = old >> 27
        return (
            (x >> rotation | x << ((32 - rotation) & 31)) & 0xFFFFFFFF
        ) / 0x100000000

    def candidate():
        while True:
            x, y = random() * 2 - 1, random() * 2 - 1
            if x * x + y * y < 1:
                z = 2 * np.sqrt(1 - x * x - y * y) * (-1 if random() < 0.5 else 1)
                vector = np.array([x, y, z])
                vector *= random() / np.linalg.norm(vector)
                vector[2] = abs(vector[2])
                return vector

    vectors = [candidate()]
    for _ in range(1, 32):
        candidates = np.array([candidate() for _ in range(10)])
        distance = np.linalg.norm(
            candidates[:, None] - np.asarray(vectors)[None], axis=2
        ).min(axis=1)
        vectors.append(candidates[np.argmax(distance)])
    scale = 0.1 + 0.9 * ((np.arange(32) + 1) / 32) ** 2
    return (np.array(vectors) * scale[:, None]).astype(np.float32)
