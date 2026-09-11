"""Pairwise metric panels with PNG export and residue selection."""

from io import BytesIO

import numpy as np
from PIL import Image

from .themes import scale


def metric(data, atoms=()):
    value = data.get("predicted_aligned_error", data.get("matrix"))
    if value is None and "values" in data:
        n = int(data.get("size", max(int(i) for i in data["values"]) + 1))
        value = np.full((n, n), np.nan)
        for i, row in data["values"].items():
            for j, v in row.items():
                value[int(i), int(j)] = float(v)
    array = np.asarray(value, float)
    if (
        array.ndim != 2
        or array.shape[0] != array.shape[1]
        or not array.size
        or np.isinf(array).any()
    ):
        raise ValueError(
            "Pairwise metric requires a square matrix; null/NaN means missing"
        )
    residues = []
    seen = set()
    for atom in atoms:
        key = (atom.model, atom.segi, atom.chain, atom.resi)
        if key not in seen:
            seen.add(key)
            residues.append(atom)
    if residues and len(residues) != len(array):
        raise ValueError("Pairwise matrix size does not match selected residue count")
    return {
        "matrix": array,
        "residues": residues,
        "title": str(data.get("title", "Pairwise metric")),
        "max": float(data.get("max_predicted_aligned_error", np.nanmax(array))),
    }


def raster(panel, size=800):
    array = panel["matrix"]
    colors = scale(array.ravel(), [0x00441B, 0xF7FCF5], [0, panel["max"]], 0xE2E2E2)
    image = Image.fromarray(np.uint8(colors.reshape(*array.shape, 3) * 255))
    return image.resize((size, size), Image.Resampling.NEAREST)


def show(panel, cmd):
    from pymol.Qt import QtGui, QtWidgets

    if QtWidgets.QApplication.instance() is None:
        return None

    class Plot(QtWidgets.QLabel):
        def mousePressEvent(self, event):
            point = event.position() if hasattr(event, "position") else event.pos()
            n = len(panel["matrix"])
            i = min(n - 1, int(point.y() / self.height() * n))
            j = min(n - 1, int(point.x() / self.width() * n))
            self.setToolTip(f"{i + 1}, {j + 1}: {panel['matrix'][i, j]:.3g}")
            residues = panel["residues"]
            if residues:
                from .source import atom_selection

                with atom_selection(
                    cmd, [(residues[k].model, residues[k].index) for k in (i, j)]
                ) as sele:
                    cmd.select("sele", f"byres ({sele})")

    dialog = QtWidgets.QDialog()
    dialog.setWindowTitle(panel["title"])
    layout = QtWidgets.QVBoxLayout(dialog)
    label = Plot()
    stream = BytesIO()
    raster(panel, 600).save(stream, format="PNG")
    pixmap = QtGui.QPixmap()
    pixmap.loadFromData(stream.getvalue())
    label.setPixmap(pixmap)
    label.setScaledContents(True)
    label.setMinimumSize(250, 250)
    layout.addWidget(label)
    dialog.resize(620, 650)
    dialog.show()
    return dialog
