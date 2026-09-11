# molstar_style_in_pymol

An independent Python implementation of Mol*-inspired visualization for PyMOL 3.1.
It provides a single `molstar_style` command, local scientific inputs, interactive
OpenGL drawing, and native ray interoperability. Node.js, Mol*, CueMol, and network
services are not runtime dependencies.

[English guide](docs/guide.md) · [日本語ガイド](docs/ja/guide.md) · [Coverage and differences](docs/coverage.md) · [Gallery](docs/gallery.md)

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

## Gallery

[All 82 examples, with GPU and ray comparison](docs/gallery.md).

These previews show the interactive GPU renderer at `quality=high`. Protein
views use crambin (PDB 1CRN); other specimens and synthetic inputs are identified
in the full gallery. Click a preview for its 1200 x 900 image.

<table>
  <tr>
    <td width="33%" align="center"><a href="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/cartoon-gpu.png"><img src="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/cartoon-gpu-thumb.png" alt="cartoon, GPU" width="260"></a><br><code>cartoon</code></td>
    <td width="33%" align="center"><a href="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/nucleic-gpu.png"><img src="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/nucleic-gpu-thumb.png" alt="nucleic, GPU" width="260"></a><br><code>nucleic</code></td>
    <td width="33%" align="center"><a href="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/ball-and-stick-gpu.png"><img src="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/ball-and-stick-gpu-thumb.png" alt="ball-and-stick, GPU" width="260"></a><br><code>ball-and-stick</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/carbohydrate-gpu.png"><img src="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/carbohydrate-gpu-thumb.png" alt="carbohydrate, GPU" width="260"></a><br><code>carbohydrate</code></td>
    <td width="33%" align="center"><a href="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/molecular-surface-gpu.png"><img src="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/molecular-surface-gpu-thumb.png" alt="molecular-surface, GPU" width="260"></a><br><code>molecular-surface</code></td>
    <td width="33%" align="center"><a href="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/direct-volume-gpu.png"><img src="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/direct-volume-gpu-thumb.png" alt="direct-volume, GPU" width="260"></a><br><code>direct-volume</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/particle-fibers-gpu.png"><img src="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/particle-fibers-gpu-thumb.png" alt="particle-fibers, GPU" width="260"></a><br><code>particle-fibers</code></td>
    <td width="33%" align="center"><a href="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/interactions-gpu.png"><img src="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/interactions-gpu-thumb.png" alt="interactions, GPU" width="260"></a><br><code>interactions</code></td>
    <td width="33%" align="center"><a href="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/orbital-gpu.png"><img src="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/orbital-gpu-thumb.png" alt="orbital, GPU" width="260"></a><br><code>orbital</code></td>
  </tr>
  <tr>
    <td width="33%" align="center"><a href="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/metallic-gpu.png"><img src="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/metallic-gpu-thumb.png" alt="metallic, GPU" width="260"></a><br><code>metallic</code></td>
    <td width="33%" align="center"><a href="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/outline-gpu.png"><img src="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/outline-gpu-thumb.png" alt="outline, GPU" width="260"></a><br><code>outline</code></td>
    <td width="33%" align="center"><a href="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/illustrative-gpu.png"><img src="https://github.com/th2ch-g/molstar_style_in_pymol/releases/download/gallery/illustrative-gpu-thumb.png" alt="illustrative, GPU" width="260"></a><br><code>illustrative</code></td>
  </tr>
</table>

Generated images are stored as release assets and remain ignored in Git.
The [rendering script](scripts/render_gallery.py) and
[regeneration instructions](docs/gallery.md#regenerate) are included.

## Install

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
