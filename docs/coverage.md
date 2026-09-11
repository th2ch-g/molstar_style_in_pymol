# Coverage and approximation contract

Reference: Mol* `5b1b54ed03b03936041f514b33b8bb129b774d37`.
The tests exercise generated geometry, ownership, finite values, reference formulas,
local formats, source restoration, real GPU output, and native ray images. A listed
visualization is a drawing capability with the documented local schema, not full
compatibility with every Mol* parameter, upstream algorithm, data service, or UI.

| Family | Implemented drawing | Differences from Mol* |
| --- | --- | --- |
| Structure (17) | cartoon, backbone, ball-and-stick, blob-surface, carbohydrate, ellipsoid, gaussian-surface, gaussian-volume, label, line, molecular-surface, orientation, plane, point, putty, spacefill, polyhedron | Mesh tessellation instead of WebGL impostors; polymer transport/interpolation and ribbon sections differ; line/point glyphs have world-space sizes |
| Cartoon visuals | trace, gaps, nucleotide block/ring/atomic fill/bonds/elements, direction wedge | Native nucleic ring unions use convex planar hulls; sugar/base geometry is approximate |
| Atomic visuals | unit/structure sphere, ellipsoid, intra/inter bonds, points/crosses | Unit identity is object/chain/segment; missing anisotropic tensors become spheres |
| Surfaces | molecular mesh/wireframe, Gaussian mesh/wireframe/volume, blob mesh/wireframe | Molecular surface uses Euclidean erosion of solvent-accessible voxels, not exact analytical rolling-probe tori; blob uses fitted group ellipsoids or regularized harmonic radial surfaces |
| Carbohydrates | SNFG symbols and internal/terminal links; 2028 recognized residue names | Polygonal approximations and face-color partitions; not exact Mol* triangulation |
| Volume (5) | direct-volume, dot, isosurface, segment, slice | Gaussian and marching-cubes grid discretization differs; dedicated ray volume uses sampled transparent planes with finite resolution; standard ray uses three pre-integrated density projections |
| Particles (4) | spacefill, orientation, fibers, target | Local particle schema and explicit positions/frames; no simulation or browser interaction manager |
| Loci/shape | distance, angle, dihedral, labels, orientation, plane, unit cell | World-space text and shapes; source atom picking and PAE pair selection are provided |
| Computed properties | geometric interactions, Shrake–Rupley ASA, explicit cross-link bounds; optional overlap clashes | Geometry-based estimates differ from Mol* feature typing and validation services; annotation inputs are preferred for external evaluations |
| Structural extensions | membrane, assembly symmetry, confal pyramids, NtC tube, clashes, tunnels | Consume local results; no ANVIL optimization, symmetry server, DNATCO assignment, or tunnel search is inferred |
| Alpha orbitals | real spherical Gaussian basis L=0..4; orbital and occupied density surfaces | CPU collocation; same reference polynomials and order conventions; finite grid/cutoff/isovalue sampling |
| Meshes/segmentation | indexed JSON/NPZ and mesh CIF/BCIF; integer label grids | Local geometry only; no segmentation editing UI or streaming server |
| Kinemage/G3D | primitive list drawing and local compressed genome coordinates | No browser controls; Kinemage ribbons use tube segments |
| MVS | local structure/component/representation/color/opacity/transform/label trees and snapshot selection | Explicit subset of MVS protocol; annotations/complex selectors require normalized schemas; unsupported nodes raise errors |
| Quality/charge views | geometry, density, RCI, pLDDT, QMEAN, partial charges | Use supplied annotations; no quality/charge model is evaluated automatically |
| Pairwise metrics | PAE/dense/sparse matrix panel, missing-data color, residue picking, PNG export | Separate Qt panel; one matrix export at a time |
| Themes | 39 built-in color names, 6 size names, extension annotation palettes | Mol* palettes/tables with independently applied interpolation/Lab adjustment; themes apply only to compatible data families |
| Materials | matte, plastic, glossy, metallic; roughness/metalness/bumpiness | Local Blinn-style shader approximation; transparent CGO and ray use baked vertex tones |
| Screen effects | occlusion, outline, shadow, cel, xray, unlit, bloom, DOF, illumination, antialias | Layer-local screen-space effects; illumination is not progressive path tracing; ray post-effects use approximate mesh-depth sampling |
| Background | gradient, radial gradient, local image, six local skybox faces | Prepared bitmap approximation; skybox updates through the maintenance timer; dedicated ray export composites backgrounds |
| Clipping | six half-spaces | Shared uncapped clipped triangles and masked volume voxels; not a full clipping-object editor |

## Validation

`tests/test_geometry.py` covers every standard visual branch and every drawing
family. `tests/test_data_themes.py` checks theme outputs, orbital/ASA/measurement
formulas, map axes/units, BinaryCIF masks, and local MVS. `tests/test_lifecycle.py`
uses real PyMOL to verify overlaps, failures, manual deletion, sessions, and states.
`tests/check_real.py --gui` records GPU and ray images and fails on empty output or
OpenGL errors. `--visuals` tests the standard subvisuals separately. Regenerable
reports, images, and benchmark results are deliberately excluded from Git.

The 500-residue, 100-state benchmark uses a repeated protein backbone with synthetic
displacements and a 1280×720 viewport. It reports preparation time, both CPU payload
sizes, bounded GPU cache use, rotation/state-switch FPS, actual movie callback state
rate, and peak RSS. Re-run on the target machine; frame rates are not universal.
