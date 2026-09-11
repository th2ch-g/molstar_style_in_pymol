"""Renderer-independent geometry shared by interactive and native-ray output."""

from dataclasses import dataclass, field

import numpy as np

from .mesh import Mesh


@dataclass
class PickTarget:
    model: str = ""
    index: int = 0
    label: str = ""


@dataclass
class Piece:
    mesh: Mesh
    atoms: tuple = ()
    edges: np.ndarray = field(default_factory=lambda: np.empty((0, 12), np.float32))
    unlit: bool = False


@dataclass
class Volume:
    grid: object
    transfer: np.ndarray
    opacity: float = 1.0
    step: float = 0.5
    pixels: object = None
    lookup: object = None

    def __post_init__(self):
        self.prepare()

    def prepare(self):
        self.pixels = np.ascontiguousarray(self.grid.values.transpose(2, 1, 0))
        lo, hi = self.transfer[0, 0], self.transfer[-1, 0]
        values = np.linspace(lo, hi if hi > lo else lo + 1, 1024)
        self.lookup = np.stack(
            [
                np.interp(values, self.transfer[:, 0], self.transfer[:, i])
                for i in range(1, 5)
            ],
            axis=-1,
        ).astype(np.float32)


@dataclass
class Geometry:
    pieces: list = field(default_factory=list)
    volumes: list = field(default_factory=list)
    panels: list = field(default_factory=list)
    metadata: dict = field(default_factory=dict)

    def add(self, mesh, atoms=(), unlit=False):
        if len(mesh.faces):
            self.pieces.append(Piece(mesh, tuple(atoms), unlit=unlit))
        return self

    def extend(self, other):
        self.pieces.extend(other.pieces)
        self.volumes.extend(other.volumes)
        self.panels.extend(other.panels)
        self.metadata.update(other.metadata)
        return self

    @property
    def nbytes(self):
        return sum(p.mesh.nbytes + p.edges.nbytes for p in self.pieces) + sum(
            v.grid.nbytes + v.transfer.nbytes for v in self.volumes
        )
