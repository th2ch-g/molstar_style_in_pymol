# Visualization gallery

[README](../README.md#gallery) · [English guide](guide.md) · [日本語ガイド](ja/guide.md) · [Coverage and differences](coverage.md)

See the [actual Mol* comparison](fidelity.md) ([日本語](ja/fidelity.md)) for
matched-camera 8GNG, material and DNA images with measured geometry and pixel errors.

82 examples cover every drawable named style, plus the nucleic-acid cartoon
alias. `empty` intentionally removes the managed geometry and has no preview.
Click an image for its 1200 x 900 PNG. Left: interactive GPU; right: dedicated
`molstar_style ray`. Both use `quality=high` and a shared camera for each example.

Protein examples use [crambin, PDB 1CRN](https://www.rcsb.org/structure/1CRN).
Atomic views use PyMOL's tryptophan fragment; DNA uses its `fnab` builder.
Glycans, coordination, volumes, particles, measurements, and extension annotations
use deterministic synthetic fixtures. Quality, validation, partial-charge, and
anisotropic-displacement values are demonstrations, not experimental results.
The ellipsoid example uses heavy atoms with a synthetic positive-definite tensor
in square angstroms. Volume dots use `sizeFactor=0.12` to separate neighboring dots;
the command default is 1. Other examples use the style defaults.
The pairwise metric is a 2D panel exported by both operations, not a ray-traced
3D object. GPU and ray differences are described in the coverage document.

Examples that require annotations or grids use the local dictionaries in the
[rendering script](../tests/render_gallery.py) and [fixtures](../tests/samples.py).
Their style names alone are not sufficient to reproduce those inputs; see the
guide for `data` schemas. Protein examples use the same orientation and framing.

Published GPU and ray PNGs are versioned in `docs/gallery/` and referenced
with relative paths, so the gallery is available directly from the checkout
and on GitHub. Images are rendered locally; no CI or release service is required.
The validation manifest stays in the ignored `.cache/gallery/` directory.

[Structure](#structure) · [Volume](#volume) · [Particles](#particles) · [Measurements and shapes](#measurements-and-shapes) · [Extensions](#extensions) · [Materials](#materials) · [Effects](#effects) · [Presets](#presets)

## Structure

| Style / specimen | GPU | Ray |
| --- | --- | --- |
| `molstar_style cartoon`<br>Crambin (PDB 1CRN) | <a href="gallery/cartoon-gpu.png"><img src="gallery/cartoon-gpu.png" alt="cartoon, GPU" width="300"></a> | <a href="gallery/cartoon-ray.png"><img src="gallery/cartoon-ray.png" alt="cartoon, RAY" width="300"></a> |
| `molstar_style backbone`<br>Crambin (PDB 1CRN) | <a href="gallery/backbone-gpu.png"><img src="gallery/backbone-gpu.png" alt="backbone, GPU" width="300"></a> | <a href="gallery/backbone-ray.png"><img src="gallery/backbone-ray.png" alt="backbone, RAY" width="300"></a> |
| `molstar_style ball-and-stick`<br>PyMOL tryptophan fragment | <a href="gallery/ball-and-stick-gpu.png"><img src="gallery/ball-and-stick-gpu.png" alt="ball-and-stick, GPU" width="300"></a> | <a href="gallery/ball-and-stick-ray.png"><img src="gallery/ball-and-stick-ray.png" alt="ball-and-stick, RAY" width="300"></a> |
| `molstar_style blob-surface`<br>Crambin (PDB 1CRN) | <a href="gallery/blob-surface-gpu.png"><img src="gallery/blob-surface-gpu.png" alt="blob-surface, GPU" width="300"></a> | <a href="gallery/blob-surface-ray.png"><img src="gallery/blob-surface-ray.png" alt="blob-surface, RAY" width="300"></a> |
| `molstar_style carbohydrate`<br>Synthetic demonstration data | <a href="gallery/carbohydrate-gpu.png"><img src="gallery/carbohydrate-gpu.png" alt="carbohydrate, GPU" width="300"></a> | <a href="gallery/carbohydrate-ray.png"><img src="gallery/carbohydrate-ray.png" alt="carbohydrate, RAY" width="300"></a> |
| `molstar_style ellipsoid`<br>PyMOL tryptophan fragment | <a href="gallery/ellipsoid-gpu.png"><img src="gallery/ellipsoid-gpu.png" alt="ellipsoid, GPU" width="300"></a> | <a href="gallery/ellipsoid-ray.png"><img src="gallery/ellipsoid-ray.png" alt="ellipsoid, RAY" width="300"></a> |
| `molstar_style gaussian-surface`<br>Crambin (PDB 1CRN) | <a href="gallery/gaussian-surface-gpu.png"><img src="gallery/gaussian-surface-gpu.png" alt="gaussian-surface, GPU" width="300"></a> | <a href="gallery/gaussian-surface-ray.png"><img src="gallery/gaussian-surface-ray.png" alt="gaussian-surface, RAY" width="300"></a> |
| `molstar_style gaussian-volume`<br>Crambin (PDB 1CRN) | <a href="gallery/gaussian-volume-gpu.png"><img src="gallery/gaussian-volume-gpu.png" alt="gaussian-volume, GPU" width="300"></a> | <a href="gallery/gaussian-volume-ray.png"><img src="gallery/gaussian-volume-ray.png" alt="gaussian-volume, RAY" width="300"></a> |
| `molstar_style label`<br>PyMOL tryptophan fragment | <a href="gallery/label-gpu.png"><img src="gallery/label-gpu.png" alt="label, GPU" width="300"></a> | <a href="gallery/label-ray.png"><img src="gallery/label-ray.png" alt="label, RAY" width="300"></a> |
| `molstar_style line`<br>PyMOL tryptophan fragment | <a href="gallery/line-gpu.png"><img src="gallery/line-gpu.png" alt="line, GPU" width="300"></a> | <a href="gallery/line-ray.png"><img src="gallery/line-ray.png" alt="line, RAY" width="300"></a> |
| `molstar_style molecular-surface`<br>Crambin (PDB 1CRN) | <a href="gallery/molecular-surface-gpu.png"><img src="gallery/molecular-surface-gpu.png" alt="molecular-surface, GPU" width="300"></a> | <a href="gallery/molecular-surface-ray.png"><img src="gallery/molecular-surface-ray.png" alt="molecular-surface, RAY" width="300"></a> |
| `molstar_style orientation`<br>Crambin (PDB 1CRN) | <a href="gallery/orientation-gpu.png"><img src="gallery/orientation-gpu.png" alt="orientation, GPU" width="300"></a> | <a href="gallery/orientation-ray.png"><img src="gallery/orientation-ray.png" alt="orientation, RAY" width="300"></a> |
| `molstar_style plane`<br>Crambin (PDB 1CRN) | <a href="gallery/plane-gpu.png"><img src="gallery/plane-gpu.png" alt="plane, GPU" width="300"></a> | <a href="gallery/plane-ray.png"><img src="gallery/plane-ray.png" alt="plane, RAY" width="300"></a> |
| `molstar_style point`<br>PyMOL tryptophan fragment | <a href="gallery/point-gpu.png"><img src="gallery/point-gpu.png" alt="point, GPU" width="300"></a> | <a href="gallery/point-ray.png"><img src="gallery/point-ray.png" alt="point, RAY" width="300"></a> |
| `molstar_style putty`<br>Crambin (PDB 1CRN) | <a href="gallery/putty-gpu.png"><img src="gallery/putty-gpu.png" alt="putty, GPU" width="300"></a> | <a href="gallery/putty-ray.png"><img src="gallery/putty-ray.png" alt="putty, RAY" width="300"></a> |
| `molstar_style spacefill`<br>PyMOL tryptophan fragment | <a href="gallery/spacefill-gpu.png"><img src="gallery/spacefill-gpu.png" alt="spacefill, GPU" width="300"></a> | <a href="gallery/spacefill-ray.png"><img src="gallery/spacefill-ray.png" alt="spacefill, RAY" width="300"></a> |
| `molstar_style polyhedron`<br>Synthetic demonstration data | <a href="gallery/polyhedron-gpu.png"><img src="gallery/polyhedron-gpu.png" alt="polyhedron, GPU" width="300"></a> | <a href="gallery/polyhedron-ray.png"><img src="gallery/polyhedron-ray.png" alt="polyhedron, RAY" width="300"></a> |
| `molstar_style nucleic`<br>DNA built with PyMOL fnab | <a href="gallery/nucleic-gpu.png"><img src="gallery/nucleic-gpu.png" alt="nucleic, GPU" width="300"></a> | <a href="gallery/nucleic-ray.png"><img src="gallery/nucleic-ray.png" alt="nucleic, RAY" width="300"></a> |

## Volume

| Style / specimen | GPU | Ray |
| --- | --- | --- |
| `molstar_style direct-volume`<br>Synthetic demonstration data | <a href="gallery/direct-volume-gpu.png"><img src="gallery/direct-volume-gpu.png" alt="direct-volume, GPU" width="300"></a> | <a href="gallery/direct-volume-ray.png"><img src="gallery/direct-volume-ray.png" alt="direct-volume, RAY" width="300"></a> |
| `molstar_style dot`<br>Synthetic demonstration data | <a href="gallery/dot-gpu.png"><img src="gallery/dot-gpu.png" alt="dot, GPU" width="300"></a> | <a href="gallery/dot-ray.png"><img src="gallery/dot-ray.png" alt="dot, RAY" width="300"></a> |
| `molstar_style isosurface`<br>Synthetic demonstration data | <a href="gallery/isosurface-gpu.png"><img src="gallery/isosurface-gpu.png" alt="isosurface, GPU" width="300"></a> | <a href="gallery/isosurface-ray.png"><img src="gallery/isosurface-ray.png" alt="isosurface, RAY" width="300"></a> |
| `molstar_style segment`<br>Synthetic demonstration data | <a href="gallery/segment-gpu.png"><img src="gallery/segment-gpu.png" alt="segment, GPU" width="300"></a> | <a href="gallery/segment-ray.png"><img src="gallery/segment-ray.png" alt="segment, RAY" width="300"></a> |
| `molstar_style slice`<br>Synthetic demonstration data | <a href="gallery/slice-gpu.png"><img src="gallery/slice-gpu.png" alt="slice, GPU" width="300"></a> | <a href="gallery/slice-ray.png"><img src="gallery/slice-ray.png" alt="slice, RAY" width="300"></a> |

## Particles

| Style / specimen | GPU | Ray |
| --- | --- | --- |
| `molstar_style particle-spacefill`<br>Synthetic demonstration data | <a href="gallery/particle-spacefill-gpu.png"><img src="gallery/particle-spacefill-gpu.png" alt="particle-spacefill, GPU" width="300"></a> | <a href="gallery/particle-spacefill-ray.png"><img src="gallery/particle-spacefill-ray.png" alt="particle-spacefill, RAY" width="300"></a> |
| `molstar_style particle-orientation`<br>Synthetic demonstration data | <a href="gallery/particle-orientation-gpu.png"><img src="gallery/particle-orientation-gpu.png" alt="particle-orientation, GPU" width="300"></a> | <a href="gallery/particle-orientation-ray.png"><img src="gallery/particle-orientation-ray.png" alt="particle-orientation, RAY" width="300"></a> |
| `molstar_style particle-fibers`<br>Synthetic demonstration data | <a href="gallery/particle-fibers-gpu.png"><img src="gallery/particle-fibers-gpu.png" alt="particle-fibers, GPU" width="300"></a> | <a href="gallery/particle-fibers-ray.png"><img src="gallery/particle-fibers-ray.png" alt="particle-fibers, RAY" width="300"></a> |
| `molstar_style particle-target`<br>Synthetic demonstration data | <a href="gallery/particle-target-gpu.png"><img src="gallery/particle-target-gpu.png" alt="particle-target, GPU" width="300"></a> | <a href="gallery/particle-target-ray.png"><img src="gallery/particle-target-ray.png" alt="particle-target, RAY" width="300"></a> |

## Measurements and shapes

| Style / specimen | GPU | Ray |
| --- | --- | --- |
| `molstar_style distance`<br>Synthetic demonstration data | <a href="gallery/distance-gpu.png"><img src="gallery/distance-gpu.png" alt="distance, GPU" width="300"></a> | <a href="gallery/distance-ray.png"><img src="gallery/distance-ray.png" alt="distance, RAY" width="300"></a> |
| `molstar_style angle`<br>Synthetic demonstration data | <a href="gallery/angle-gpu.png"><img src="gallery/angle-gpu.png" alt="angle, GPU" width="300"></a> | <a href="gallery/angle-ray.png"><img src="gallery/angle-ray.png" alt="angle, RAY" width="300"></a> |
| `molstar_style dihedral`<br>Synthetic demonstration data | <a href="gallery/dihedral-gpu.png"><img src="gallery/dihedral-gpu.png" alt="dihedral, GPU" width="300"></a> | <a href="gallery/dihedral-ray.png"><img src="gallery/dihedral-ray.png" alt="dihedral, RAY" width="300"></a> |
| `molstar_style shape-label`<br>Synthetic demonstration data | <a href="gallery/shape-label-gpu.png"><img src="gallery/shape-label-gpu.png" alt="shape-label, GPU" width="300"></a> | <a href="gallery/shape-label-ray.png"><img src="gallery/shape-label-ray.png" alt="shape-label, RAY" width="300"></a> |
| `molstar_style shape-orientation`<br>Synthetic demonstration data | <a href="gallery/shape-orientation-gpu.png"><img src="gallery/shape-orientation-gpu.png" alt="shape-orientation, GPU" width="300"></a> | <a href="gallery/shape-orientation-ray.png"><img src="gallery/shape-orientation-ray.png" alt="shape-orientation, RAY" width="300"></a> |
| `molstar_style shape-plane`<br>Synthetic demonstration data | <a href="gallery/shape-plane-gpu.png"><img src="gallery/shape-plane-gpu.png" alt="shape-plane, GPU" width="300"></a> | <a href="gallery/shape-plane-ray.png"><img src="gallery/shape-plane-ray.png" alt="shape-plane, RAY" width="300"></a> |
| `molstar_style unitcell`<br>Synthetic demonstration data | <a href="gallery/unitcell-gpu.png"><img src="gallery/unitcell-gpu.png" alt="unitcell, GPU" width="300"></a> | <a href="gallery/unitcell-ray.png"><img src="gallery/unitcell-ray.png" alt="unitcell, RAY" width="300"></a> |

## Extensions

| Style / specimen | GPU | Ray |
| --- | --- | --- |
| `molstar_style interactions`<br>Synthetic demonstration data | <a href="gallery/interactions-gpu.png"><img src="gallery/interactions-gpu.png" alt="interactions, GPU" width="300"></a> | <a href="gallery/interactions-ray.png"><img src="gallery/interactions-ray.png" alt="interactions, RAY" width="300"></a> |
| `molstar_style cross-link-restraint`<br>Synthetic demonstration data | <a href="gallery/cross-link-restraint-gpu.png"><img src="gallery/cross-link-restraint-gpu.png" alt="cross-link-restraint, GPU" width="300"></a> | <a href="gallery/cross-link-restraint-ray.png"><img src="gallery/cross-link-restraint-ray.png" alt="cross-link-restraint, RAY" width="300"></a> |
| `molstar_style membrane-orientation`<br>Synthetic demonstration data | <a href="gallery/membrane-orientation-gpu.png"><img src="gallery/membrane-orientation-gpu.png" alt="membrane-orientation, GPU" width="300"></a> | <a href="gallery/membrane-orientation-ray.png"><img src="gallery/membrane-orientation-ray.png" alt="membrane-orientation, RAY" width="300"></a> |
| `molstar_style assembly-symmetry`<br>Synthetic demonstration data | <a href="gallery/assembly-symmetry-gpu.png"><img src="gallery/assembly-symmetry-gpu.png" alt="assembly-symmetry, GPU" width="300"></a> | <a href="gallery/assembly-symmetry-ray.png"><img src="gallery/assembly-symmetry-ray.png" alt="assembly-symmetry, RAY" width="300"></a> |
| `molstar_style confal-pyramids`<br>Synthetic demonstration data | <a href="gallery/confal-pyramids-gpu.png"><img src="gallery/confal-pyramids-gpu.png" alt="confal-pyramids, GPU" width="300"></a> | <a href="gallery/confal-pyramids-ray.png"><img src="gallery/confal-pyramids-ray.png" alt="confal-pyramids, RAY" width="300"></a> |
| `molstar_style ntc-tube`<br>Synthetic demonstration data | <a href="gallery/ntc-tube-gpu.png"><img src="gallery/ntc-tube-gpu.png" alt="ntc-tube, GPU" width="300"></a> | <a href="gallery/ntc-tube-ray.png"><img src="gallery/ntc-tube-ray.png" alt="ntc-tube, RAY" width="300"></a> |
| `molstar_style clashes`<br>Synthetic demonstration data | <a href="gallery/clashes-gpu.png"><img src="gallery/clashes-gpu.png" alt="clashes, GPU" width="300"></a> | <a href="gallery/clashes-ray.png"><img src="gallery/clashes-ray.png" alt="clashes, RAY" width="300"></a> |
| `molstar_style orbital`<br>Synthetic demonstration data | <a href="gallery/orbital-gpu.png"><img src="gallery/orbital-gpu.png" alt="orbital, GPU" width="300"></a> | <a href="gallery/orbital-ray.png"><img src="gallery/orbital-ray.png" alt="orbital, RAY" width="300"></a> |
| `molstar_style orbital-density`<br>Synthetic demonstration data | <a href="gallery/orbital-density-gpu.png"><img src="gallery/orbital-density-gpu.png" alt="orbital-density, GPU" width="300"></a> | <a href="gallery/orbital-density-ray.png"><img src="gallery/orbital-density-ray.png" alt="orbital-density, RAY" width="300"></a> |
| `molstar_style tunnel`<br>Synthetic demonstration data | <a href="gallery/tunnel-gpu.png"><img src="gallery/tunnel-gpu.png" alt="tunnel, GPU" width="300"></a> | <a href="gallery/tunnel-ray.png"><img src="gallery/tunnel-ray.png" alt="tunnel, RAY" width="300"></a> |
| `molstar_style mesh`<br>Synthetic demonstration data | <a href="gallery/mesh-gpu.png"><img src="gallery/mesh-gpu.png" alt="mesh, GPU" width="300"></a> | <a href="gallery/mesh-ray.png"><img src="gallery/mesh-ray.png" alt="mesh, RAY" width="300"></a> |
| `molstar_style kinemage`<br>Synthetic demonstration data | <a href="gallery/kinemage-gpu.png"><img src="gallery/kinemage-gpu.png" alt="kinemage, GPU" width="300"></a> | <a href="gallery/kinemage-ray.png"><img src="gallery/kinemage-ray.png" alt="kinemage, RAY" width="300"></a> |
| `molstar_style g3d`<br>Synthetic demonstration data | <a href="gallery/g3d-gpu.png"><img src="gallery/g3d-gpu.png" alt="g3d, GPU" width="300"></a> | <a href="gallery/g3d-ray.png"><img src="gallery/g3d-ray.png" alt="g3d, RAY" width="300"></a> |
| `molstar_style mvs`<br>Crambin (PDB 1CRN) | <a href="gallery/mvs-gpu.png"><img src="gallery/mvs-gpu.png" alt="mvs, GPU" width="300"></a> | <a href="gallery/mvs-ray.png"><img src="gallery/mvs-ray.png" alt="mvs, RAY" width="300"></a> |
| `molstar_style pairwise-metric`<br>Synthetic demonstration data | <a href="gallery/pairwise-metric-gpu.png"><img src="gallery/pairwise-metric-gpu.png" alt="pairwise-metric, GPU" width="300"></a> | <a href="gallery/pairwise-metric-ray.png"><img src="gallery/pairwise-metric-ray.png" alt="pairwise-metric, RAY" width="300"></a> |
| `molstar_style annotation-label`<br>Synthetic demonstration data | <a href="gallery/annotation-label-gpu.png"><img src="gallery/annotation-label-gpu.png" alt="annotation-label, GPU" width="300"></a> | <a href="gallery/annotation-label-ray.png"><img src="gallery/annotation-label-ray.png" alt="annotation-label, RAY" width="300"></a> |
| `molstar_style custom-label`<br>Synthetic demonstration data | <a href="gallery/custom-label-gpu.png"><img src="gallery/custom-label-gpu.png" alt="custom-label, GPU" width="300"></a> | <a href="gallery/custom-label-ray.png"><img src="gallery/custom-label-ray.png" alt="custom-label, RAY" width="300"></a> |

## Materials

| Style / specimen | GPU | Ray |
| --- | --- | --- |
| `molstar_style matte`<br>Crambin (PDB 1CRN) | <a href="gallery/matte-gpu.png"><img src="gallery/matte-gpu.png" alt="matte, GPU" width="300"></a> | <a href="gallery/matte-ray.png"><img src="gallery/matte-ray.png" alt="matte, RAY" width="300"></a> |
| `molstar_style plastic`<br>Crambin (PDB 1CRN) | <a href="gallery/plastic-gpu.png"><img src="gallery/plastic-gpu.png" alt="plastic, GPU" width="300"></a> | <a href="gallery/plastic-ray.png"><img src="gallery/plastic-ray.png" alt="plastic, RAY" width="300"></a> |
| `molstar_style glossy`<br>Crambin (PDB 1CRN) | <a href="gallery/glossy-gpu.png"><img src="gallery/glossy-gpu.png" alt="glossy, GPU" width="300"></a> | <a href="gallery/glossy-ray.png"><img src="gallery/glossy-ray.png" alt="glossy, RAY" width="300"></a> |
| `molstar_style metallic`<br>Crambin (PDB 1CRN) | <a href="gallery/metallic-gpu.png"><img src="gallery/metallic-gpu.png" alt="metallic, GPU" width="300"></a> | <a href="gallery/metallic-ray.png"><img src="gallery/metallic-ray.png" alt="metallic, RAY" width="300"></a> |

## Effects

| Style / specimen | GPU | Ray |
| --- | --- | --- |
| `molstar_style outline`<br>Crambin (PDB 1CRN) | <a href="gallery/outline-gpu.png"><img src="gallery/outline-gpu.png" alt="outline, GPU" width="300"></a> | <a href="gallery/outline-ray.png"><img src="gallery/outline-ray.png" alt="outline, RAY" width="300"></a> |
| `molstar_style occlusion`<br>Crambin (PDB 1CRN) | <a href="gallery/occlusion-gpu.png"><img src="gallery/occlusion-gpu.png" alt="occlusion, GPU" width="300"></a> | <a href="gallery/occlusion-ray.png"><img src="gallery/occlusion-ray.png" alt="occlusion, RAY" width="300"></a> |
| `molstar_style shadow`<br>Crambin (PDB 1CRN) | <a href="gallery/shadow-gpu.png"><img src="gallery/shadow-gpu.png" alt="shadow, GPU" width="300"></a> | <a href="gallery/shadow-ray.png"><img src="gallery/shadow-ray.png" alt="shadow, RAY" width="300"></a> |
| `molstar_style cel`<br>Crambin (PDB 1CRN) | <a href="gallery/cel-gpu.png"><img src="gallery/cel-gpu.png" alt="cel, GPU" width="300"></a> | <a href="gallery/cel-ray.png"><img src="gallery/cel-ray.png" alt="cel, RAY" width="300"></a> |
| `molstar_style xray`<br>Crambin (PDB 1CRN) | <a href="gallery/xray-gpu.png"><img src="gallery/xray-gpu.png" alt="xray, GPU" width="300"></a> | <a href="gallery/xray-ray.png"><img src="gallery/xray-ray.png" alt="xray, RAY" width="300"></a> |
| `molstar_style unlit`<br>Crambin (PDB 1CRN) | <a href="gallery/unlit-gpu.png"><img src="gallery/unlit-gpu.png" alt="unlit, GPU" width="300"></a> | <a href="gallery/unlit-ray.png"><img src="gallery/unlit-ray.png" alt="unlit, RAY" width="300"></a> |
| `molstar_style bloom`<br>Crambin (PDB 1CRN) | <a href="gallery/bloom-gpu.png"><img src="gallery/bloom-gpu.png" alt="bloom, GPU" width="300"></a> | <a href="gallery/bloom-ray.png"><img src="gallery/bloom-ray.png" alt="bloom, RAY" width="300"></a> |
| `molstar_style dof`<br>Crambin (PDB 1CRN) | <a href="gallery/dof-gpu.png"><img src="gallery/dof-gpu.png" alt="dof, GPU" width="300"></a> | <a href="gallery/dof-ray.png"><img src="gallery/dof-ray.png" alt="dof, RAY" width="300"></a> |
| `molstar_style illumination`<br>Crambin (PDB 1CRN) | <a href="gallery/illumination-gpu.png"><img src="gallery/illumination-gpu.png" alt="illumination, GPU" width="300"></a> | <a href="gallery/illumination-ray.png"><img src="gallery/illumination-ray.png" alt="illumination, RAY" width="300"></a> |
| `molstar_style background`<br>Crambin (PDB 1CRN) | <a href="gallery/background-gpu.png"><img src="gallery/background-gpu.png" alt="background, GPU" width="300"></a> | <a href="gallery/background-ray.png"><img src="gallery/background-ray.png" alt="background, RAY" width="300"></a> |
| `molstar_style antialias`<br>Crambin (PDB 1CRN) | <a href="gallery/antialias-gpu.png"><img src="gallery/antialias-gpu.png" alt="antialias, GPU" width="300"></a> | <a href="gallery/antialias-ray.png"><img src="gallery/antialias-ray.png" alt="antialias, RAY" width="300"></a> |

## Presets

| Style / specimen | GPU | Ray |
| --- | --- | --- |
| `molstar_style default`<br>Crambin (PDB 1CRN) | <a href="gallery/default-gpu.png"><img src="gallery/default-gpu.png" alt="default, GPU" width="300"></a> | <a href="gallery/default-ray.png"><img src="gallery/default-ray.png" alt="default, RAY" width="300"></a> |
| `molstar_style auto`<br>Crambin (PDB 1CRN) | <a href="gallery/auto-gpu.png"><img src="gallery/auto-gpu.png" alt="auto, GPU" width="300"></a> | <a href="gallery/auto-ray.png"><img src="gallery/auto-ray.png" alt="auto, RAY" width="300"></a> |
| `molstar_style polymer-and-ligand`<br>Crambin (PDB 1CRN) | <a href="gallery/polymer-and-ligand-gpu.png"><img src="gallery/polymer-and-ligand-gpu.png" alt="polymer-and-ligand, GPU" width="300"></a> | <a href="gallery/polymer-and-ligand-ray.png"><img src="gallery/polymer-and-ligand-ray.png" alt="polymer-and-ligand, RAY" width="300"></a> |
| `molstar_style protein-and-nucleic`<br>Crambin (PDB 1CRN) | <a href="gallery/protein-and-nucleic-gpu.png"><img src="gallery/protein-and-nucleic-gpu.png" alt="protein-and-nucleic, GPU" width="300"></a> | <a href="gallery/protein-and-nucleic-ray.png"><img src="gallery/protein-and-nucleic-ray.png" alt="protein-and-nucleic, RAY" width="300"></a> |
| `molstar_style polymer-cartoon`<br>Crambin (PDB 1CRN) | <a href="gallery/polymer-cartoon-gpu.png"><img src="gallery/polymer-cartoon-gpu.png" alt="polymer-cartoon, GPU" width="300"></a> | <a href="gallery/polymer-cartoon-ray.png"><img src="gallery/polymer-cartoon-ray.png" alt="polymer-cartoon, RAY" width="300"></a> |
| `molstar_style atomic-detail`<br>Crambin (PDB 1CRN) | <a href="gallery/atomic-detail-gpu.png"><img src="gallery/atomic-detail-gpu.png" alt="atomic-detail, GPU" width="300"></a> | <a href="gallery/atomic-detail-ray.png"><img src="gallery/atomic-detail-ray.png" alt="atomic-detail, RAY" width="300"></a> |
| `molstar_style coarse-surface`<br>Crambin (PDB 1CRN) | <a href="gallery/coarse-surface-gpu.png"><img src="gallery/coarse-surface-gpu.png" alt="coarse-surface, GPU" width="300"></a> | <a href="gallery/coarse-surface-ray.png"><img src="gallery/coarse-surface-ray.png" alt="coarse-surface, RAY" width="300"></a> |
| `molstar_style illustrative`<br>Crambin (PDB 1CRN) | <a href="gallery/illustrative-gpu.png"><img src="gallery/illustrative-gpu.png" alt="illustrative, GPU" width="300"></a> | <a href="gallery/illustrative-ray.png"><img src="gallery/illustrative-ray.png" alt="illustrative, RAY" width="300"></a> |
| `molstar_style auto-lod`<br>Crambin (PDB 1CRN) | <a href="gallery/auto-lod-gpu.png"><img src="gallery/auto-lod-gpu.png" alt="auto-lod, GPU" width="300"></a> | <a href="gallery/auto-lod-ray.png"><img src="gallery/auto-lod-ray.png" alt="auto-lod, RAY" width="300"></a> |
| `molstar_style mesoscale`<br>Crambin (PDB 1CRN) | <a href="gallery/mesoscale-gpu.png"><img src="gallery/mesoscale-gpu.png" alt="mesoscale, GPU" width="300"></a> | <a href="gallery/mesoscale-ray.png"><img src="gallery/mesoscale-ray.png" alt="mesoscale, RAY" width="300"></a> |
| `molstar_style validation-geometry`<br>1CRN with synthetic annotation values | <a href="gallery/validation-geometry-gpu.png"><img src="gallery/validation-geometry-gpu.png" alt="validation-geometry, GPU" width="300"></a> | <a href="gallery/validation-geometry-ray.png"><img src="gallery/validation-geometry-ray.png" alt="validation-geometry, RAY" width="300"></a> |
| `molstar_style validation-density`<br>1CRN with synthetic annotation values | <a href="gallery/validation-density-gpu.png"><img src="gallery/validation-density-gpu.png" alt="validation-density, GPU" width="300"></a> | <a href="gallery/validation-density-ray.png"><img src="gallery/validation-density-ray.png" alt="validation-density, RAY" width="300"></a> |
| `molstar_style validation-rci`<br>1CRN with synthetic annotation values | <a href="gallery/validation-rci-gpu.png"><img src="gallery/validation-rci-gpu.png" alt="validation-rci, GPU" width="300"></a> | <a href="gallery/validation-rci-ray.png"><img src="gallery/validation-rci-ray.png" alt="validation-rci, RAY" width="300"></a> |
| `molstar_style quality-plddt`<br>1CRN with synthetic annotation values | <a href="gallery/quality-plddt-gpu.png"><img src="gallery/quality-plddt-gpu.png" alt="quality-plddt, GPU" width="300"></a> | <a href="gallery/quality-plddt-ray.png"><img src="gallery/quality-plddt-ray.png" alt="quality-plddt, RAY" width="300"></a> |
| `molstar_style quality-qmean`<br>1CRN with synthetic annotation values | <a href="gallery/quality-qmean-gpu.png"><img src="gallery/quality-qmean-gpu.png" alt="quality-qmean, GPU" width="300"></a> | <a href="gallery/quality-qmean-ray.png"><img src="gallery/quality-qmean-ray.png" alt="quality-qmean, RAY" width="300"></a> |
| `molstar_style partial-charges`<br>1CRN with synthetic annotation values | <a href="gallery/partial-charges-gpu.png"><img src="gallery/partial-charges-gpu.png" alt="partial-charges, GPU" width="300"></a> | <a href="gallery/partial-charges-ray.png"><img src="gallery/partial-charges-ray.png" alt="partial-charges, RAY" width="300"></a> |

## Regenerate

The renderer accepts local inputs and performs no downloads. Run it in the
repository's PyMOL Qt environment; it opens a separate window.

```sh
pixi install --locked
uv run --no-project --python .pixi/envs/default/bin/python python \
    tests/render_gallery.py --structure path/to/1crn.cif
```

Published images go to `docs/gallery/`. Use `--only cartoon nucleic`
for selected examples or `--output .cache/gallery-preview` for a temporary run.
The ignored `.cache/gallery/manifest.json` records the source revision, input
checksum, image checksums, and nonempty-image validation. Use `--manifest` to
choose another report location. On Linux without a desktop, prefix the command
with `xvfb-run -a`.
