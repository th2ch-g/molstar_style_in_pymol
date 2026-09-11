# Mol*-inspired visualization in PyMOL

[日本語](ja/guide.md) · [Coverage](coverage.md) · [GPU/ray gallery](gallery.md)

`molstar_style style, selection=all` applies a named managed view. The default
is `polymer-and-ligand`: polymer cartoon, ligand/ion sticks and spheres, water
with partial opacity, and SNFG glycans. `protein-and-nucleic` isolates polymers.
The package independently implements the visualization families in Mol* revision
`5b1b54ed03b03936041f514b33b8bb129b774d37`, including extension geometry.

## Command

```text
molstar_style [style], selection=all, representation=auto, color=auto, quality=medium, name=molstar
molstar_style list
molstar_style help
molstar_style refresh, name=all
molstar_style reset, name=all
molstar_style png, filename=figure.png, width=1600, height=1200
molstar_style ray, filename=figure_ray.png, width=1600, height=1200
```

`representation` overrides geometry when applying a material or effect. `color=keep`
uses the original PyMOL atom colors. `list` prints all style, color, and size names.
`quality` accepts `lowest`, `lower`, `low`, `medium`, `high`, `higher`, `highest`,
`auto`, and `custom`. Use `params` for custom resolution or tessellation.
`transparency=0..1` multiplies each layer's opacity; `keep` retains preset opacity.
`state=0` prepares all loaded states; a positive state number creates one fixed view.

`data` accepts a local file path or Python dictionary. `params` accepts a local
JSON path or Python dictionary. PyMOL's comma parser is best used with JSON files;
use Python calls for dictionaries. Inputs and referenced assets must be local.
No command downloads structures, annotations, or assets. Missing required scientific
annotations cause an error before replacing the current managed view.

```python
from molstar_style_in_pymol import molstar_style
molstar_style('cartoon', 'chain A', params={
    'tubularHelices': True, 'aspectRatio': 5,
    'visuals': ['polymer-trace', 'polymer-gap', 'nucleotide-ring'],
    'material': {'metalness': 0.2, 'roughness': 0.4, 'bumpiness': 0},
    'postprocessing': {'occlusion': True, 'outline': True},
})
```

## Representation parameters

| Family | Parameters |
| --- | --- |
| Atoms | `sizeFactor`, `ignoreHydrogens`, `sizeAspectRatio`, `multipleBonds` (`off`, `symmetric`, `offset`), `aromaticBonds`, `visuals` |
| Ellipsoids | Anisotropic atom tensors; `probability` defaults to 0.5; isotropic atoms use spheres |
| Cartoon/backbone/putty | `sizeFactor`, `aspectRatio`, `arrowFactor`, `tubularHelices`, `linearSegments`, `radialSegments`, `bfactorScale`, `visuals` |
| Surface | `resolution` in angstroms, `probeRadius`, `radiusOffset`, `smoothness`, `isoValue`, mesh/wireframe `visuals` |
| Blob | `blobSize`, `method` (`grid`, `clustering`), `clusterIterations`, `shape` (`ellipsoid`, `spherical-harmonics`), `degree`, `regularization` |
| Volume isosurface/dot | `isoValue` (absolute number or `{"kind":"relative","relativeValue":1}`), `visuals`, `showWireframe`, `sizeFactor`; dots: `stride`, `maxPoints` |
| Slice | `dimension` (`x`, `y`, `z`), `index` or `relativeIndex`, `colorList`, `domain` |
| Direct volume | `transferFunction`: increasing `[value, color, opacity]` rows, or normalized `controlPoints`: `[fraction, opacity]`; `step` in angstroms |
| Segments | Integer label grid, `segments` list, `smoothness`, `colorList` |
| Labels | `level` (`chain`, `residue`, `element`), `sizeFactor`; custom labels: `text`, `textSize` |
| Measurements | `sizeFactor`, `textSize`, `label`, `arcScale`; positions must define nondegenerate geometry |
| Clipping | `clipPlanes`: up to six `[nx,ny,nz,offset]` half-spaces; keeps `n dot position + offset >= 0` |

See `reference.json` (package data) for every standard `visuals` name. Both
unit and structure visual variants use the same independently generated mesh;
inter-unit bonds use distinct object/chain/segment identities. A visual requiring
absent geometry, such as polymer gaps in a continuous chain, can be empty.

## Color, size, and annotations

All 39 built-in color names and six size names are exposed, together with extension
color names. Molecular themes use `colorParams`; size themes use `sizeTheme` and
`sizeParams`. Intrinsic themes read atom identities, element, chain, secondary
structure, occupancy, B-factor, charges, or van der Waals radii from PyMOL. The
physical source coordinates and properties are never recolored or rewritten.

Annotations use one of these explicit alignments:

```json
{"residues": [{"chain": "A", "resi": "10", "plddt": 94.5}]}
```

```json
{"atom_data": {"entity_id": ["1", "1"], "partial_charge": [-0.3, 0.3]}}
```

`atom_data` arrays must match the selected atom order from `cmd.get_model`.
Residue keys include optional `model` and `segi`; duplicate matches are rejected.
Missing residue values use the missing-data color. Scalar extension fields are
`plddt`, `qmean`, `geometry_quality`, `density_fit`, `random_coil_index`,
`pdbe_structure_quality_report`, and `sb_ncbr_partial_charges`. Categorical fields
use the color name with hyphens replaced by underscores. Entity/operator themes
require actual entity/operator metadata. Local mmCIF/BCIF model-archive QA tables
and entity-source tables are normalized automatically; additional validation
formats can be converted to the residue schema above.

`accessible-surface-area` calculates Shrake–Rupley area in square angstroms.
`external-volume` samples a supplied local grid at atom coordinates.
`external-structure` uses explicit `external_coordinates` and `external_colors`.
Particle colors use particle entity/compartment/hierarchy/index/attribute fields;
volume colors use value/instance/segment fields. Inapplicable themes fail explicitly.

## Local input schemas

Files: CCP4/MRC, Gaussian Cube, OpenDX, NumPy NPZ (`values`, `transform`),
CIF/BCIF density-server grids, JSON, MVS JSON, binary G3D, and Kinemage text.
PyMOL map objects can be passed as `selection`. Grid transforms map voxel indices
to angstrom coordinates, including nonorthogonal cells and axis permutations.
NPZ loading disables pickle. `data` dictionaries can contain a `Grid` directly.

| Display | Required JSON fields |
| --- | --- |
| Measurements | `positions`: two/three/four XYZ rows; alternatively `indices` into selected atoms |
| Custom/annotation labels | `positions`, `text`; selected atom center is the default position |
| Mesh | `vertices`, integer triangle `faces`; optional `colors`, `color`, `opacity`, `transform`, `label`; multiple items via `meshes` |
| Mesh BCIF | `mesh`, `mesh_vertex`, `mesh_triangle` tables |
| Particles | `particles`: rows with `position`, `radius`, optional `quaternion` in XYZW order, `axes`, `scale`, `entity`, `compartment`, `hierarchy`, `color`, `label` |
| Fibers/targets | Particle rows also require `points` / `target` |
| Cross-links | `cross_links`: rows with `indices: [i,j]`, `lower`, `upper`; colors report distance violations |
| Annotated contacts/clashes | `interactions` / `clashes`: rows with `indices`, `type`, optional `positions`, `color` |
| Membrane | `membrane`: `center`, `normal`, `thickness`, `radius` |
| Assembly symmetry | `symmetry.axes`: `start`, `end`, `order`; optional `symmetry.cage.vertices` and `edges` |
| DNATCO | `steps`: `class`, `score`, `positions`; confal order is O3-prime, P, OP1, OP2, O5-prime |
| Tunnel | `positions`, positive `radii`; multiple tunnels via `tunnels` |
| Orbital | Mol* `basis.atoms` with Bohr centers and shells (`exponents`, `angularMomentum`, `coefficients`); `orbitals` with `alpha`, `occupancy`, `energy` |
| Pairwise metric/PAE | Square `predicted_aligned_error` or `matrix`, optional title and maximum; sparse nested `values` is also accepted |
| Unit cell | `cell: [a,b,c,alpha,beta,gamma]`, optional origin; can read PyMOL symmetry |

`interactions` without annotations computes geometric hydrogen bonds, optional weak
hydrogen bonds, hydrophobic/ionic/halogen/metal contacts, pi stacking, and cation–pi
contacts. Explicit hydrogens constrain donor angles. These are geometric estimates,
not a protonation or energy calculation. `clashes` requires annotations unless
`params={"compute": true}` explicitly requests van der Waals overlap estimates.
Membrane, symmetry, DNATCO, and validation displays consume supplied results;
they do not infer those scientific annotations from geometry.

Orbital evaluation follows Mol* real solid harmonics L=0..4 and its `gaussian`,
`cca`, and `cca-reverse` coefficient orders. Density sums occupied squared orbitals.
A precomputed Cube grid also works. Kinemage handles vectors, ribbons, triangles,
balls, spheres, dots, labels, and words. G3D reads local compressed resolution
blocks, with haplotype, chromosome, and region filters.

MVS supports local structure trees (`root/download/parse/structure/component/`
`representation/color/opacity/transform/label`) and selecting a snapshot from
multiple-state documents. URLs must be local paths. Complex annotation selectors
must first be normalized to the explicit local annotation/mesh schemas above;
unsupported nodes fail instead of silently being omitted. MVS is a visualization
input, not a browser, server, animation editor, or full Mol* application replacement.

PAE opens a Qt panel; clicking a cell selects its two residues when matrix rows
match the selected residues. PNG/ray operations export a matrix image when the
managed view consists only of a panel. Headless matrix export is supported.

## Effects and output

Materials: `matte`, `plastic`, `glossy`, `metallic`, or a dictionary with
`metalness`, `roughness`, `bumpiness` in `[0,1]`. Effects: `occlusion`, `outline`,
`shadow`, `cel`, `xray`, `unlit`, `bloom`, `dof`, `illumination`, `antialias`.
`postprocessing` values can be booleans or dictionaries with `strength`/`intensity`.
`focus` is normalized depth for depth-of-field effects. Default occlusion is scoped
to the managed layer. Illumination approximates indirect lighting in screen space;
this is not Mol*'s progressive path tracer. Ray effects use mesh-depth samples.

`background` accepts Mol*'s `variant` dictionary with `horizontalGradient`,
`radialGradient`, `image`, or `skybox`. Images and six cube faces must be local.
Rotation, blur, saturation, and lightness are prepared outside drawing callbacks.
The existing background remains unchanged unless a background is explicitly requested.

Opaque geometry uses GLSL/VBOs. Transparent mesh bodies use native CGO with baked
vertex lighting. Volumes use actual scalar-field ray marching in OpenGL. Standard
PyMOL `ray` and `png, ray=1` include retained meshes and three pre-integrated density projections;
these standard commands are not patched. Default PyMOL transparency mode 2
does not accumulate overlapping transparent surfaces, so the retained projections
provide a coarse preview without changing global settings. Use `molstar_style ray` for camera-aligned
density sampling, explicit mesh outlines, background compositing, and image effects.
Different renderers and sampling produce visible differences; see the coverage table.

## Lifecycle and performance

Names can overlap on the same atoms; reset restores an atom's original representation
when the final owner is removed. Reapplying a name prepares and loads the replacement
before releasing the previous view. Failed preparation/loading preserves the old view.
Deleting a managed group or its generated object triggers cleanup through the Qt
maintenance timer or the next command. Sessions store restoration records; reloading
rebuilds CPU/GPU geometry, and missing local files leave native sources restored.
Use `refresh` after modifying source coordinates, topology, colors, or annotations.

All source states are prepared outside OpenGL callbacks. `cache_mb=2048` bounds
geometry plus native float CGO payloads; allocator/Python overhead and source snapshots
are additional. `gpu_cache_mb=256` bounds VBO and volume-texture caches; framebuffers
are additional viewport-dependent allocations. The command reports preparation time
and geometry size. The reproducible benchmark reports native payload, GPU cache, and
peak process RSS; its coordinates are repeated protein backbones with synthetic motion,
not an MD simulation. GPU, platform, selection, and quality affect playback speed.

```sh
pixi run test
pixi run check
pixi run uv run --no-project python tests/check_real.py --gui
pixi run uv run --no-project python tests/check_real.py --gui --visuals
pixi run uv run --no-project python tests/check_real.py --gui --benchmark --structure local_protein.cif
```

GUI checks create a separate PyMOL process. They never reuse a live user session.
Gallery images and machine-specific benchmark reports belong in `.cache/`.

ASCII labels use the bundled vector font. Unicode labels use the Qt font system
when a GUI is present; headless Unicode rendering requires `params.font` pointing
to a local TTF/OTF font containing the requested glyphs. Text is tessellated once
and is shared by GPU and ray output.
