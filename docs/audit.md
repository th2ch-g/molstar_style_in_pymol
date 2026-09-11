# Geometry audit

[English guide](guide.md) · [日本語](ja/guide.md) · [Gallery](gallery.md) · [Coverage](coverage.md)

Reference: Mol* revision `5b1b54ed03b03936041f514b33b8bb129b774d37`.
This review compares the local implementation with that source, then checks
physical invariants numerically and renders each family in PyMOL. Rendering a
nonempty image alone does not establish agreement with Mol*.

## Corrections

| Area | Incorrect behavior | Corrected behavior and regression evidence |
| --- | --- | --- |
| Ellipsoid bonds | Physical radii and aspect ratio 2/3 made carbon bonds 1.1333 Å thick | Uniform size 1 and aspect ratio 0.1 give 0.1 Å; tests also vary theme, aspect ratio, and factor |
| Ellipsoid atoms | Missing tensors became van der Waals spheres; reflected eigenvectors reversed winding | Missing tensors are omitted; missing selections fail transactionally; tensor equation, analytic normals, isotropic limit, reflected winding, PDB ANISOU units/order are checked |
| Ellipsoid scale | Exact 50% chi-square quantile was incorrectly presented as the Mol* default | Default uses Mol*'s literal 1.5958 multiplier and absolute eigenvalues; explicit `probability` remains a separate exact-quantile option |
| Backbone, putty | Backbone radius 0.2; putty used a square-root B-factor formula; polymer size themes were ignored | Backbone radius 0.3; putty radius `0.2 * (0.2 + 0.1 * B)`; straight-chain bounds and explicit size themes are checked |
| Cartoon subvisuals | Default nucleic atomic bonds/elements/fill omitted; nucleic connectors leaked into polymer-only visuals | Include default nucleic components; draw base connectors only when a nucleotide visual is requested |
| Lines and multiple bonds | Boolean attenuation used as a radius; lone atoms disappeared; crosses replaced points; arbitrary perpendicular offsets could leave a ring plane | World-space glyph conversion is explicit; separate point/cross handling; lone crosses; reference mesh radius/spacing formulas with neighboring-atom offset planes |
| Carbohydrate visuals | Internal/terminal links were not separated; reflected symbol frames were possible | Filter each link by its endpoints, preserve right-handed symbol frames, honor explicit molecular colors |
| Surfaces and orbitals | Triangle winding opposed normals; negative orbital lobes had inward normals | Correct winding after marching cubes, grid reflections, and negative isovalues; tests compare face/vertex normals and outward directions in anisotropic grids |
| Structure plane | Plain uniformly colored fitted rectangle | Atom-colored sphere cross sections, physical radii, offset/rotation/plane controls, cutout; an analytic off-center sphere section checks its radius |
| Orientation | Structure default drew axes, box, and ellipsoid together; size factor ignored; PCA box centered at atom mean | Structure defaults to ellipsoid; shape defaults to box; use projected min/max midpoint and scale around it; rotation/translation invariants are checked |
| Volume colors | Several accepted color options were ignored, including molecular Gaussian-volume colors | Uniform colors affect all five volume families; slices modulate scalar intensity; Gaussian color grids feed GPU, dedicated ray, and retained projections |
| Volume dots | Negative thresholds selected the wrong side; radius default 0.12 | Negative thresholds select lower values; default radius 1; gallery explicitly requests 0.12 |
| Angle/dihedral | Fixed-radius outline without default sector; valid 0/180-degree angles rejected | Relative arm scaling, sectors and visual controls; sector area and collinear angle limits are checked; undefined dihedrals still fail |
| Particles | Explicit scale discarded radius; orientation used radius-dependent arrows; fibers were joined straight cylinders; target meant an arrow endpoint | Multiplicative scale, fixed-length orientation axes, interpolated fiber tube, centered and rotated local shape/structure/volume instances |
| Confal pyramids | Phosphate was used as the dividing vertex; arbitrary category colors | O3'/O5' midpoint and upper/lower NtC reference colors; vertex location and AB01 colors are checked |
| Mesh CIF/BCIF | Vertex IDs assumed to be contiguous array offsets | Resolve vertex IDs per mesh and reject duplicates/missing references; reordered sparse IDs are checked |
| Native ray and effects | CGO triangle opcodes used the same vertex order; back-facing planes lost their color; depth lookup could index outside the image; vertex-only depth estimates caused speckled occlusion | Correct opcode-specific winding, temporarily enable two-sided lighting, rasterize triangle interiors within image bounds, decode PyMOL background color IDs, and preserve untouched flat-background pixels. Real PyMOL tests render a colored plane from both sides and check setting restoration |

`tests/test_reference_geometry.py` contains the numerical regressions;
`tests/test_lifecycle.py` checks real PyMOL tensor ingestion and failure recovery.
`tests/check_real.py --gui --visuals` checks each standard subvisual separately.
`tests/render_gallery.py` renders all 82 gallery entries locally in both renderers.
The gallery, including its synthetic inputs, is not an experimental validation dataset.

Local results for this correction: 253 passing pytest cases; 55 standard subvisuals
rendered in GPU/ray; 82 gallery pairs (164 PNGs, 1200 × 900); Ruff and source/wheel
builds passed. No CI rendering is used. These checks establish the behaviors above,
not complete upstream parameter or image equivalence.

## Family review and remaining approximations

| Families reviewed | Remaining differences / scope |
| --- | --- |
| cartoon, backbone, putty, nucleic alias | Independent Hermite interpolation, frame transport, gap detection, ribbon profiles, and convex base-ring geometry; not all Mol* control points or profile parameters |
| ball-and-stick, line, point, spacefill, ellipsoid | Meshes replace impostors; line/point widths are in Å rather than pixels; bond reference selection, trimming, aromatic perception, symmetry/unit identity and filtering differ. Nonpositive ellipsoid eigenvalues follow Mol* absolute-value handling; this does not validate the scientific quality of ADPs |
| molecular-surface, gaussian-surface, gaussian-volume, blob-surface | Gaussian grid evaluation and marching-cubes sampling are independent; molecular surface is a voxel erosion approximation, not exact rolling-probe geometry; blob grouping/harmonic fits differ; Gaussian colors use the nearest atom |
| carbohydrate, polyhedron, orientation, plane, label | SNFG tables are retained, but shape tessellation/face partitions and ring detection differ; coordination uses a distance cutoff; PCA degeneracies, plane rasterization and trimming, and text metrics are independent |
| direct-volume, dot, isosurface, segment, slice | Local scalar grids and transforms; no periodic wrapping, flood-fill editor, GPU impostors, or oblique volume slices. Transfer-function opacity, finite GPU step counts, and sampled native-ray density differ |
| particle-spacefill, particle-orientation, particle-fibers, particle-target | Explicit local particle/target schemas; shape target center defaults to the bounding-box midpoint; structure blob targets use a Gaussian envelope; no simulation, streaming, per-target dynamic LOD or Mol* instancing backend |
| distance, angle, dihedral, shape-label, shape-orientation, shape-plane, unitcell | Numerical measurements and cell validity are checked; labels use fixed world-space geometry, not screen-facing billboards; line thickness, orientation axis radii and fitted-plane boundaries differ |
| interactions, cross-link-restraint, clashes | Explicit annotation geometry and bound violations; computed contacts are geometric estimates rather than Mol* chemical feature typing. Clash glyphs use sphere/line markers rather than Mol* disk glyphs |
| membrane-orientation, assembly-symmetry | Draw supplied planes, normals, symmetry axes and cage edges; no ANVIL optimization or external symmetry assignment; labels and glyph sizes differ |
| confal-pyramids, ntc-tube | Reference pyramid coordinates/half colors; NtC tube consumes a supplied path and class colors, without Mol* residue/step boundary construction or DNATCO assignment |
| orbital, orbital-density | Reference real solid harmonics L=0..4 and coefficient orders; density is the occupied sum of squared orbitals. CPU grid/cutoff sampling and default 15%-of-maximum isovalue differ from cumulative-probability threshold selection |
| tunnel, mesh, kinemage, g3d | Local profiles, meshes and genomic coordinates; tunnel interpolation and kinemage ribbons are independent; no tunnel search, streaming or browser controls |
| mvs, annotation-label, custom-label, pairwise-metric | Supported local MVS subset rejects unsupported nodes; annotation labels require normalized inputs; PAE is a Qt matrix panel with residue-pair selection, not a Mol* browser panel |
| default, auto, empty, polymer-and-ligand, protein-and-nucleic, polymer-cartoon, atomic-detail, coarse-surface, illustrative, auto-lod, mesoscale | Independent preset composition; auto/auto-lod/mesoscale use static atom-count thresholds, without Mol* structure-size categories or distance-dependent LOD |
| validation-geometry, validation-density, validation-rci, quality-plddt, quality-qmean, partial-charges | Color supplied aligned annotations; do not calculate validation metrics or partial charges |
| matte, plastic, glossy, metallic | Reference material parameter values feed a local shader; no physically equivalent Mol* BRDF; transparent CGO/ray lighting uses baked tones |
| outline, occlusion, shadow, cel, xray, unlit, bloom, dof, illumination, background, antialias | Independent screen-space/ray effects; illumination is not a path tracer. Backgrounds use local bitmap preparation; sampling and compositing cannot guarantee pixel equivalence |

All named families remain drawable through the documented interface. This is a
Mol*-inspired PyMOL implementation, not a complete parameter-compatible port.
Differences above are explicit limitations, not evidence of exact Mol* agreement.

## Primary source locations

- [Ellipsoid representation defaults](https://github.com/molstar/molstar/blob/5b1b54ed03b03936041f514b33b8bb129b774d37/src/mol-repr/structure/representation/ellipsoid.ts), [tensor geometry](https://github.com/molstar/molstar/blob/5b1b54ed03b03936041f514b33b8bb129b774d37/src/mol-repr/structure/visual/ellipsoid-mesh.ts).
- [Uncertainty size theme](https://github.com/molstar/molstar/blob/5b1b54ed03b03936041f514b33b8bb129b774d37/src/mol-theme/size/uncertainty.ts), [bond geometry](https://github.com/molstar/molstar/blob/5b1b54ed03b03936041f514b33b8bb129b774d37/src/mol-repr/structure/visual/util/link.ts).
- [Structure plane](https://github.com/molstar/molstar/blob/5b1b54ed03b03936041f514b33b8bb129b774d37/src/mol-repr/structure/visual/plane-image.ts), [volume slice](https://github.com/molstar/molstar/blob/5b1b54ed03b03936041f514b33b8bb129b774d37/src/mol-repr/volume/slice.ts).
- [Particle target transforms](https://github.com/molstar/molstar/blob/5b1b54ed03b03936041f514b33b8bb129b774d37/src/mol-repr/particles/representation/target/common.ts), [confal pyramids](https://github.com/molstar/molstar/blob/5b1b54ed03b03936041f514b33b8bb129b774d37/src/extensions/dnatco/confal-pyramids/representation.ts).
