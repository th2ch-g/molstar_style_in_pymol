"""Compatibility-profile GLSL renderer; never calls the PyMOL command API."""

from collections import OrderedDict
from ctypes import c_void_p
from dataclasses import dataclass, field
from importlib.resources import files

import numpy as np


@dataclass
class Drawing:
    pieces: list
    profile: object
    edge_color: tuple
    pool: object
    volumes: list = field(default_factory=list)
    name: str = ""
    state: int = 1
    matrices: object = None
    error: str = ""
    draws: int = 0
    anchors: object = None
    keys: tuple = ()
    selected: object = None
    extent: object = field(default_factory=lambda: [[0, 0, 0], [0, 0, 0]])

    def __post_init__(self):
        arrays = [p.mesh.vertices for p in self.pieces if len(p.mesh.vertices)]
        if arrays:
            self.extent = [
                np.min([a.min(axis=0) for a in arrays], axis=0).tolist(),
                np.max([a.max(axis=0) for a in arrays], axis=0).tolist(),
            ]

    def get_extent(self):
        return self.extent

    def __getstate__(self):
        # Session restore tasks rebuild geometry; GL handles are context-local.
        return {"extent": self.extent}

    def __setstate__(self, state):
        self.extent = state["extent"]
        self.pool = None

    def __call__(self):
        if self.pool is None or self.error:
            return
        try:
            self.pool.draw(self)
            self.draws += 1
        except Exception as exc:
            self.error = str(exc)
            print(
                f" molstar_style: OpenGL drawing failed: {exc}. Use refresh after correcting the OpenGL context."
            )


class Pool:
    """Bounded VBO cache shared by all states in a managed view."""

    def __init__(self, budget_mb=256):
        self.budget = int(budget_mb * 1024**2)
        self.buffers = OrderedDict()
        self.programs = {}
        self.context = None
        self.bytes = 0
        self.show_selection = True
        self.framebuffer = None
        self.volume_buffers = OrderedDict()

    def clear(self):
        from OpenGL import GL as gl
        from OpenGL import contextdata

        current = contextdata.getContext()
        if current == self.context:
            for handles, _, _ in self.buffers.values():
                gl.glDeleteBuffers(len(handles), handles)
            for program in self.programs.values():
                gl.glDeleteProgram(program)
            for handles, _ in self.volume_buffers.values():
                gl.glDeleteTextures(handles)
            if self.framebuffer is not None:
                self.framebuffer.close()
        self.framebuffer = None
        self.volume_buffers.clear()
        self.buffers.clear()
        self.programs.clear()
        self.bytes = 0
        self.context = current

    def program(self, name):
        from OpenGL import GL as gl
        from OpenGL.GL.shaders import compileProgram, compileShader

        if name not in self.programs:
            root = files(__package__).joinpath("shaders")
            self.programs[name] = compileProgram(
                compileShader(
                    root.joinpath(name + ".vert").read_text(), gl.GL_VERTEX_SHADER
                ),
                compileShader(
                    root.joinpath(name + ".frag").read_text(), gl.GL_FRAGMENT_SHADER
                ),
                validate=False,
            )
        return self.programs[name]

    def buffer(self, piece, edges=False):
        from OpenGL import GL as gl

        key = (id(piece.mesh), edges)
        if key in self.buffers:
            self.buffers.move_to_end(key)
            return self.buffers[key]
        if edges:
            e = piece.edges
            data = np.empty(
                (len(e), 4),
                dtype=[
                    ("position", "f4", 3),
                    ("normal_a", "i2", 3),
                    ("normal_b", "i2", 3),
                    ("opposite_side", "f4", 4),
                ],
            )
            data["position"][:, :2] = e[:, None, :3]
            data["position"][:, 2:] = e[:, None, 3:6]
            data["opposite_side"][:, :2, :3] = e[:, None, 3:6]
            data["opposite_side"][:, 2:, :3] = e[:, None, :3]
            data["opposite_side"][:, :, 3] = [-1.0, 1.0, -1.0, 1.0]
            # Signed 16-bit unit normals reduce state uploads without changing
            # positions or edge width; the vertex shader normalizes both.
            data["normal_a"] = np.rint(32767 * e[:, None, 6:9]).astype(np.int16)
            data["normal_b"] = np.rint(32767 * e[:, None, 9:12]).astype(np.int16)
            arrays = [(gl.GL_ARRAY_BUFFER, data)]
            count = data.size
        else:
            mesh = piece.mesh
            data = np.c_[
                mesh.vertices, mesh.normals, mesh.colors, mesh.alphas * mesh.opacity
            ].astype(np.float32)
            arrays = [
                (gl.GL_ARRAY_BUFFER, data),
                (gl.GL_ELEMENT_ARRAY_BUFFER, mesh.faces),
            ]
            count = mesh.faces.size
        size = sum(a.nbytes for _, a in arrays)
        self.reserve(size)
        if size > self.budget:
            raise RuntimeError(
                "A GPU buffer exceeds 256 MiB; lower quality or select fewer atoms"
            )
        while self.buffers and self.bytes + size > self.budget:
            _, (handles, _, old_size) = self.buffers.popitem(last=False)
            gl.glDeleteBuffers(len(handles), handles)
            self.bytes -= old_size
        handles = [int(gl.glGenBuffers(1)) for _ in arrays]
        try:
            for handle, (target, array) in zip(handles, arrays):
                gl.glBindBuffer(target, handle)
                gl.glBufferData(target, array.nbytes, array, gl.GL_STATIC_DRAW)
        except Exception:
            gl.glDeleteBuffers(len(handles), handles)
            raise
        value = handles, count, size
        self.buffers[key] = value
        self.bytes += size
        return value

    def reserve(self, size):
        from OpenGL import GL as gl

        if size > self.budget:
            raise RuntimeError("A GPU resource exceeds gpu_cache_mb")
        while self.bytes + size > self.budget:
            if self.buffers:
                _, (handles, _, old_size) = self.buffers.popitem(last=False)
                gl.glDeleteBuffers(len(handles), handles)
            elif self.volume_buffers:
                _, (handles, old_size) = self.volume_buffers.popitem(last=False)
                gl.glDeleteTextures(handles)
            else:
                break
            self.bytes -= old_size

    def draw(self, drawing):
        from OpenGL import GL as gl
        from OpenGL import contextdata

        from .postprocessing import Framebuffer, effect_value

        current = contextdata.getContext()
        if current != self.context:
            self.clear()
        old_program = int(gl.glGetIntegerv(gl.GL_CURRENT_PROGRAM))
        old_array = int(gl.glGetIntegerv(gl.GL_ARRAY_BUFFER_BINDING))
        old_element = int(gl.glGetIntegerv(gl.GL_ELEMENT_ARRAY_BUFFER_BINDING))
        old_client = int(gl.glGetIntegerv(gl.GL_CLIENT_ACTIVE_TEXTURE))
        old_active = int(gl.glGetIntegerv(gl.GL_ACTIVE_TEXTURE))
        modelview = np.asarray(gl.glGetDoublev(gl.GL_MODELVIEW_MATRIX)).T.copy()
        projection = np.asarray(gl.glGetDoublev(gl.GL_PROJECTION_MATRIX)).T.copy()
        viewport = np.asarray(gl.glGetIntegerv(gl.GL_VIEWPORT)).copy()
        drawing.matrices = (modelview, projection, viewport)
        gl.glPushAttrib(gl.GL_ALL_ATTRIB_BITS)
        gl.glPushClientAttrib(gl.GL_CLIENT_ALL_ATTRIB_BITS)
        framebuffer = None
        try:
            if drawing.profile.background:
                from .background import draw

                draw(drawing.profile.background, self, projection, modelview, viewport)
            if any(effect_value(v) > 0 for v in drawing.profile.effects.values()):
                if self.framebuffer is None or (
                    self.framebuffer.width,
                    self.framebuffer.height,
                ) != (int(viewport[2]), int(viewport[3])):
                    if self.framebuffer is not None:
                        self.framebuffer.close()
                    self.framebuffer = Framebuffer(int(viewport[2]), int(viewport[3]))
                framebuffer = self.framebuffer
                framebuffer.begin()
            gl.glDisable(gl.GL_LIGHTING)
            gl.glDisable(gl.GL_FOG)
            gl.glDisable(gl.GL_BLEND)
            gl.glDisable(gl.GL_CULL_FACE)
            gl.glEnable(gl.GL_DEPTH_TEST)
            gl.glDepthFunc(gl.GL_LEQUAL)
            gl.glDepthMask(True)
            gl.glEnableClientState(gl.GL_VERTEX_ARRAY)
            gl.glEnableClientState(gl.GL_NORMAL_ARRAY)
            gl.glEnableClientState(gl.GL_COLOR_ARRAY)
            gl.glClientActiveTexture(gl.GL_TEXTURE0)
            gl.glDisableClientState(gl.GL_TEXTURE_COORD_ARRAY)
            program = self.program("body")
            gl.glUseProgram(program)

            def integer(key, value):
                gl.glUniform1i(gl.glGetUniformLocation(program, key), int(value))

            material = drawing.profile.material
            gl.glUniform3f(
                gl.glGetUniformLocation(program, "material"),
                *(float(material[k]) for k in ("metalness", "roughness", "bumpiness")),
            )
            integer("cel", drawing.profile.params.get("cel", False))
            integer("xray", drawing.profile.params.get("xray", False))
            planes = np.asarray(
                drawing.profile.params.get("clipPlanes", []), np.float32
            ).reshape(-1, 4)
            integer("clipCount", len(planes))
            if len(planes):
                gl.glUniform4fv(
                    gl.glGetUniformLocation(program, "clipPlanes"), len(planes), planes
                )
            for piece in drawing.pieces:
                if (
                    piece.mesh.opacity < 0.999999
                    or np.min(piece.mesh.alphas, initial=1) < 0.999999
                ):
                    continue
                integer("unlit", piece.unlit or drawing.profile.ignore_light)
                handles, count, _ = self.buffer(piece)
                gl.glBindBuffer(gl.GL_ARRAY_BUFFER, handles[0])
                gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, handles[1])
                gl.glVertexPointer(3, gl.GL_FLOAT, 40, c_void_p(0))
                gl.glNormalPointer(gl.GL_FLOAT, 40, c_void_p(12))
                gl.glColorPointer(4, gl.GL_FLOAT, 40, c_void_p(24))
                gl.glDrawElements(
                    gl.GL_TRIANGLES, count, gl.GL_UNSIGNED_INT, c_void_p(0)
                )
            if drawing.volumes:
                from .volume_gpu import draw_volumes

                draw_volumes(self, drawing)
            if framebuffer is not None:
                framebuffer.composite(self.program("post"), drawing.profile.effects)
            if (
                self.show_selection
                and drawing.selected is not None
                and len(drawing.selected)
            ):
                gl.glUseProgram(0)
                gl.glDisable(gl.GL_DEPTH_TEST)
                gl.glBindBuffer(gl.GL_ARRAY_BUFFER, 0)
                gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, 0)
                gl.glDisableClientState(gl.GL_NORMAL_ARRAY)
                gl.glDisableClientState(gl.GL_COLOR_ARRAY)
                gl.glVertexPointer(3, gl.GL_FLOAT, 0, drawing.selected)
                gl.glColor3f(1, 0.1, 0.8)
                gl.glPointSize(7)
                gl.glDrawArrays(gl.GL_POINTS, 0, len(drawing.selected))
        finally:
            if framebuffer is not None:
                framebuffer.restore()
            gl.glPopClientAttrib()
            gl.glPopAttrib()
            gl.glUseProgram(old_program)
            gl.glBindBuffer(gl.GL_ARRAY_BUFFER, old_array)
            gl.glBindBuffer(gl.GL_ELEMENT_ARRAY_BUFFER, old_element)
            gl.glClientActiveTexture(old_client)
            gl.glActiveTexture(old_active)
