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
    color_grid: object = None

    def __post_init__(self):
        self.prepare()

    def prepare(self):
        self.pixels = np.ascontiguousarray(self.grid.values.transpose(2, 1, 0))
        if self.color_grid is not None:
            self.color_grid = np.asarray(self.color_grid, np.float32)
            if (
                self.color_grid.shape != (*self.grid.values.shape, 3)
                or not np.isfinite(self.color_grid).all()
            ):
                raise ValueError("Volume colors must match the scalar grid")
            self.pixels = np.ascontiguousarray(
                np.concatenate(
                    [self.color_grid, self.grid.values[..., None]], axis=-1
                ).transpose(2, 1, 0, 3)
            )
        lo, hi = self.transfer[0, 0], self.transfer[-1, 0]
        values = np.linspace(lo, hi if hi > lo else lo + 1, 1024)
        self.lookup = np.stack(
            [
                np.interp(values, self.transfer[:, 0], self.transfer[:, i])
                for i in range(1, 5)
            ],
            axis=-1,
        ).astype(np.float32)

    def sample_colors(self, positions):
        from scipy.ndimage import map_coordinates

        inv = np.linalg.inv(self.grid.transform)
        indices = np.asarray(positions) @ inv[:3, :3].T + inv[:3, 3]
        return np.stack(
            [
                map_coordinates(
                    self.color_grid[..., k], indices.T, order=1, mode="nearest"
                )
                for k in range(3)
            ],
            axis=-1,
        )


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
            v.grid.nbytes
            + v.transfer.nbytes
            + v.pixels.nbytes
            + (0 if v.color_grid is None else v.color_grid.nbytes)
            for v in self.volumes
        )
