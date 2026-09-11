# molstar_style_in_pymol

An independent Python implementation of Mol*-inspired visualization for PyMOL 3.1.
It provides a single `molstar_style` command, local scientific inputs, interactive
OpenGL drawing, and native ray interoperability. Node.js, Mol*, CueMol, and network
services are not runtime dependencies.

[English guide](docs/guide.md) · [日本語ガイド](docs/ja/guide.md) · [Coverage and differences](docs/coverage.md)

```text
molstar_style
molstar_style cartoon, color=secondary-structure
molstar_style glossy, representation=ball-and-stick
molstar_style isosurface, selection=density
molstar_style direct-volume, data=density.mrc
molstar_style interactions, selection=protein or ligand, name=contacts
molstar_style ray, filename=figure.png, width=1600, height=1200
molstar_style reset, name=all
```

Install into the Python environment used by PyMOL:

```sh
uv pip install --python <pymol-python> git+https://github.com/th2ch-g/molstar_style_in_pymol.git
```

Register in PyMOL's Python interpreter or startup file:

```python
from molstar_style_in_pymol import __init_plugin__
__init_plugin__()
```

The `mdtbx` `pymol_plugins` package performs this registration automatically.
Importing this package alone does not import or start PyMOL.

For development, install [pixi](https://pixi.sh) and uv, then run `pixi install`
and `pixi run test`. The pinned development interpreter is Python 3.10 with
PyMOL 3.1. Interactive rendering requires compatibility OpenGL 2.1 / GLSL 1.20;
headless PyMOL supports ray export. Generated images and caches are ignored.

The reference revision is Mol* `5b1b54ed03b03936041f514b33b8bb129b774d37`.
See [NOTICE](NOTICE) for attribution and the independent renderer's scope.
