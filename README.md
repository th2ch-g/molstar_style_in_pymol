# molstar_style_in_pymol

An independent Python implementation of Mol*-inspired visualization for PyMOL 3.1.
It provides a single `molstar_style` command, local scientific inputs, interactive
OpenGL drawing, and native ray interoperability. Node.js, Mol*, CueMol, and network
services are not runtime dependencies.

[English guide](docs/guide.md) · [日本語ガイド](docs/ja/guide.md) · [Mol* comparison](docs/fidelity.md) · [Coverage and differences](docs/coverage.md) · [Geometry audit](docs/audit.md) · [Gallery](docs/gallery.md)

Ribbon curves, cross sections, GGX materials, lighting and ambient occlusion follow
the pinned Mol* source. The [matched-camera comparison](docs/fidelity.md) runs actual
Mol* and PyMOL on 8GNG, 1CRN and DNA, with numerical geometry and image metrics.

| Mol* reference (8GNG) | PyMOL GPU | PyMOL dedicated ray |
| --- | --- | --- |
| ![Molstar 8GNG](docs/gallery/reference-molstar-8gng-color-gpu.png) | ![PyMOL GPU 8GNG](docs/gallery/reference-pymol-8gng-color-gpu.png) | ![PyMOL ray 8GNG](docs/gallery/reference-pymol-8gng-color-ray.png) |

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

All 82 examples, rendered locally at 1200 x 900 with `quality=high`. Protein
views share the same camera and use [crambin, PDB 1CRN](https://www.rcsb.org/structure/1CRN).
DNA uses PyMOL's `fnab` builder; atomic views use its tryptophan fragment.
Volumes, particles, measurements, and extension inputs are synthetic examples.

These PNGs show the interactive GPU renderer. The
[full gallery](docs/gallery.md) compares every example with native ray output
and identifies its input data. `empty` intentionally has no preview.

### Structure

<table>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/cartoon-gpu.png" alt="cartoon" width="240"><br><code>cartoon</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/backbone-gpu.png" alt="backbone" width="240"><br><code>backbone</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/ball-and-stick-gpu.png" alt="ball-and-stick" width="240"><br><code>ball-and-stick</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/blob-surface-gpu.png" alt="blob-surface" width="240"><br><code>blob-surface</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/carbohydrate-gpu.png" alt="carbohydrate" width="240"><br><code>carbohydrate</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/ellipsoid-gpu.png" alt="ellipsoid" width="240"><br><code>ellipsoid</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/gaussian-surface-gpu.png" alt="gaussian-surface" width="240"><br><code>gaussian-surface</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/gaussian-volume-gpu.png" alt="gaussian-volume" width="240"><br><code>gaussian-volume</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/label-gpu.png" alt="label" width="240"><br><code>label</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/line-gpu.png" alt="line" width="240"><br><code>line</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/molecular-surface-gpu.png" alt="molecular-surface" width="240"><br><code>molecular-surface</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/orientation-gpu.png" alt="orientation" width="240"><br><code>orientation</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/plane-gpu.png" alt="plane" width="240"><br><code>plane</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/point-gpu.png" alt="point" width="240"><br><code>point</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/putty-gpu.png" alt="putty" width="240"><br><code>putty</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/spacefill-gpu.png" alt="spacefill" width="240"><br><code>spacefill</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/polyhedron-gpu.png" alt="polyhedron" width="240"><br><code>polyhedron</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/nucleic-gpu.png" alt="nucleic" width="240"><br><code>nucleic</code></td>
  </tr>
</table>

### Volume

<table>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/direct-volume-gpu.png" alt="direct-volume" width="240"><br><code>direct-volume</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/dot-gpu.png" alt="dot" width="240"><br><code>dot</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/isosurface-gpu.png" alt="isosurface" width="240"><br><code>isosurface</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/segment-gpu.png" alt="segment" width="240"><br><code>segment</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/slice-gpu.png" alt="slice" width="240"><br><code>slice</code></td>
    <td width="33%"></td>
  </tr>
</table>

### Particles

<table>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/particle-spacefill-gpu.png" alt="particle-spacefill" width="240"><br><code>particle-spacefill</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/particle-orientation-gpu.png" alt="particle-orientation" width="240"><br><code>particle-orientation</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/particle-fibers-gpu.png" alt="particle-fibers" width="240"><br><code>particle-fibers</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/particle-target-gpu.png" alt="particle-target" width="240"><br><code>particle-target</code></td>
    <td width="33%"></td>
    <td width="33%"></td>
  </tr>
</table>

### Measurements and shapes

<table>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/distance-gpu.png" alt="distance" width="240"><br><code>distance</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/angle-gpu.png" alt="angle" width="240"><br><code>angle</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/dihedral-gpu.png" alt="dihedral" width="240"><br><code>dihedral</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/shape-label-gpu.png" alt="shape-label" width="240"><br><code>shape-label</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/shape-orientation-gpu.png" alt="shape-orientation" width="240"><br><code>shape-orientation</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/shape-plane-gpu.png" alt="shape-plane" width="240"><br><code>shape-plane</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/unitcell-gpu.png" alt="unitcell" width="240"><br><code>unitcell</code></td>
    <td width="33%"></td>
    <td width="33%"></td>
  </tr>
</table>

### Extensions

<table>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/interactions-gpu.png" alt="interactions" width="240"><br><code>interactions</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/cross-link-restraint-gpu.png" alt="cross-link-restraint" width="240"><br><code>cross-link-restraint</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/membrane-orientation-gpu.png" alt="membrane-orientation" width="240"><br><code>membrane-orientation</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/assembly-symmetry-gpu.png" alt="assembly-symmetry" width="240"><br><code>assembly-symmetry</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/confal-pyramids-gpu.png" alt="confal-pyramids" width="240"><br><code>confal-pyramids</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/ntc-tube-gpu.png" alt="ntc-tube" width="240"><br><code>ntc-tube</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/clashes-gpu.png" alt="clashes" width="240"><br><code>clashes</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/orbital-gpu.png" alt="orbital" width="240"><br><code>orbital</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/orbital-density-gpu.png" alt="orbital-density" width="240"><br><code>orbital-density</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/tunnel-gpu.png" alt="tunnel" width="240"><br><code>tunnel</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/mesh-gpu.png" alt="mesh" width="240"><br><code>mesh</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/kinemage-gpu.png" alt="kinemage" width="240"><br><code>kinemage</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/g3d-gpu.png" alt="g3d" width="240"><br><code>g3d</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/mvs-gpu.png" alt="mvs" width="240"><br><code>mvs</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/pairwise-metric-gpu.png" alt="pairwise-metric" width="240"><br><code>pairwise-metric</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/annotation-label-gpu.png" alt="annotation-label" width="240"><br><code>annotation-label</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/custom-label-gpu.png" alt="custom-label" width="240"><br><code>custom-label</code></td>
    <td width="33%"></td>
  </tr>
</table>

### Materials

<table>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/matte-gpu.png" alt="matte" width="240"><br><code>matte</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/plastic-gpu.png" alt="plastic" width="240"><br><code>plastic</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/glossy-gpu.png" alt="glossy" width="240"><br><code>glossy</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/metallic-gpu.png" alt="metallic" width="240"><br><code>metallic</code></td>
    <td width="33%"></td>
    <td width="33%"></td>
  </tr>
</table>

### Effects

<table>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/outline-gpu.png" alt="outline" width="240"><br><code>outline</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/occlusion-gpu.png" alt="occlusion" width="240"><br><code>occlusion</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/shadow-gpu.png" alt="shadow" width="240"><br><code>shadow</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/cel-gpu.png" alt="cel" width="240"><br><code>cel</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/xray-gpu.png" alt="xray" width="240"><br><code>xray</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/unlit-gpu.png" alt="unlit" width="240"><br><code>unlit</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/bloom-gpu.png" alt="bloom" width="240"><br><code>bloom</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/dof-gpu.png" alt="dof" width="240"><br><code>dof</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/illumination-gpu.png" alt="illumination" width="240"><br><code>illumination</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/background-gpu.png" alt="background" width="240"><br><code>background</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/antialias-gpu.png" alt="antialias" width="240"><br><code>antialias</code></td>
    <td width="33%"></td>
  </tr>
</table>

### Presets

<table>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/default-gpu.png" alt="default" width="240"><br><code>default</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/auto-gpu.png" alt="auto" width="240"><br><code>auto</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/polymer-and-ligand-gpu.png" alt="polymer-and-ligand" width="240"><br><code>polymer-and-ligand</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/protein-and-nucleic-gpu.png" alt="protein-and-nucleic" width="240"><br><code>protein-and-nucleic</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/polymer-cartoon-gpu.png" alt="polymer-cartoon" width="240"><br><code>polymer-cartoon</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/atomic-detail-gpu.png" alt="atomic-detail" width="240"><br><code>atomic-detail</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/coarse-surface-gpu.png" alt="coarse-surface" width="240"><br><code>coarse-surface</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/illustrative-gpu.png" alt="illustrative" width="240"><br><code>illustrative</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/auto-lod-gpu.png" alt="auto-lod" width="240"><br><code>auto-lod</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/mesoscale-gpu.png" alt="mesoscale" width="240"><br><code>mesoscale</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/validation-geometry-gpu.png" alt="validation-geometry" width="240"><br><code>validation-geometry</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/validation-density-gpu.png" alt="validation-density" width="240"><br><code>validation-density</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/validation-rci-gpu.png" alt="validation-rci" width="240"><br><code>validation-rci</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/quality-plddt-gpu.png" alt="quality-plddt" width="240"><br><code>quality-plddt</code></td>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/quality-qmean-gpu.png" alt="quality-qmean" width="240"><br><code>quality-qmean</code></td>
  </tr>
  <tr>
    <td width="33%" align="center" valign="top"><img src="docs/gallery/partial-charges-gpu.png" alt="partial-charges" width="240"><br><code>partial-charges</code></td>
    <td width="33%"></td>
    <td width="33%"></td>
  </tr>
</table>

Run `molstar_style <style>` for a named style. Examples with grids or annotations
require the local `data` inputs documented in the [full gallery](docs/gallery.md).

Regenerate the gallery in the repository's PyMOL Qt environment:

```sh
uv run --no-project --python .pixi/envs/default/bin/python python \
    tests/render_gallery.py --structure path/to/1crn.cif
```

The published PNGs are versioned for README display. Downloaded coordinates,
render manifests, and temporary validation output stay in the ignored `.cache`
directory.

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
headless PyMOL supports ray export. Published gallery PNGs are versioned;
temporary renders and caches are ignored.

The reference revision is Mol* `5b1b54ed03b03936041f514b33b8bb129b774d37`.
See [NOTICE](NOTICE) for attribution and the independent renderer's scope.
