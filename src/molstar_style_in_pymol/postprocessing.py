"""Scoped screen effects with explicit Qt framebuffer restoration."""


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

    def composite(self, program, effects):
        from OpenGL import GL as gl

        self.bind(gl.GL_FRAMEBUFFER, self.previous)
        gl.glViewport(*self.viewport)
        gl.glUseProgram(program)
        for unit, (name, texture) in enumerate(zip(("image", "depth"), self.textures)):
            gl.glActiveTexture(gl.GL_TEXTURE0 + unit)
            gl.glBindTexture(gl.GL_TEXTURE_2D, texture)
            gl.glUniform1i(gl.glGetUniformLocation(program, name), unit)
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
        gl.glBlendFunc(gl.GL_SRC_ALPHA, gl.GL_ONE_MINUS_SRC_ALPHA)
        gl.glEnable(gl.GL_DEPTH_TEST)
        gl.glDepthMask(True)
        gl.glBegin(gl.GL_QUADS)
        for x, y in ((0, 0), (1, 0), (1, 1), (0, 1)):
            gl.glTexCoord2f(x, y)
            gl.glVertex3f(x * 2 - 1, y * 2 - 1, 0)
        gl.glEnd()

    def close(self):
        from OpenGL import GL as gl

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
