# Coverage and approximation contract

[Package README](../README.md) | [Documentation](README.md) | [English guide](pymol_molstar.md) | [日本語](ja/pymol_molstar.md)

Reference: Mol* `5b1b54ed03b03936041f514b33b8bb129b774d37`.
See the [geometry audit](audit.md) and [actual Mol* comparison](fidelity.md) for
corrected formulas, measured images and remaining limitations.
The tests exercise generated geometry, ownership, finite values, reference formulas,
local formats, source restoration, real GPU output, and native ray images. A listed
visualization is a drawing capability with the documented local schema, not full
compatibility with every Mol* parameter, upstream algorithm, data service, or UI.

| Family | Implemented drawing | Differences from Mol* |
| --- | --- | --- |
| Structure (17) | cartoon, backbone, ball-and-stick, blob-surface, carbohydrate, ellipsoid, gaussian-surface, gaussian-volume, label, line, molecular-surface, orientation, plane, point, putty, spacefill, polyhedron | Mesh tessellation instead of WebGL impostors; line/point glyphs have world-space sizes |
| Cartoon visuals | Mol* residue-local trace curves and frames, elliptical/rounded/square profiles, beta arrows, helixorient tubes, gaps, nucleotide visuals, direction wedge | PyMOL H/S/coil assignments merge upstream turn categories; internal caps can differ. Cyclic/coarse polymers, gaps and base-ring hulls remain independent approximations |
| Atomic visuals | unit/structure sphere, ellipsoid, intra/inter bonds, points/crosses | Unit identity is object/chain/segment; missing anisotropic tensors are omitted, and entirely missing tensors raise an error; bond-reference selection and end trimming differ |
| Surfaces | molecular mesh/wireframe, Gaussian mesh/wireframe/volume, blob mesh/wireframe | Molecular surface uses Euclidean erosion of solvent-accessible voxels, not exact analytical rolling-probe tori; blob uses fitted group ellipsoids or regularized harmonic radial surfaces |
| Carbohydrates | SNFG symbols and internal/terminal links; 2028 recognized residue names | Polygonal approximations and face-color partitions; not exact Mol* triangulation |
| Volume (5) | direct-volume, dot, isosurface, segment, slice; pixel-wise dedicated ray density clipped by managed opaque meshes | Gaussian/grid discretization differs; unmanaged or transparent mesh scenes use sampled native planes; standard ray uses three pre-integrated projections |
| Particles (4) | spacefill, orientation axes, interpolated fibers, instanced shape/structure/volume targets | Local schema and explicit positions/frames; target center defaults to a box midpoint; structure blob targets use Gaussian surfaces; no dynamic LOD, streaming, or simulation |
| Loci/shape | distance, angle, dihedral, labels, orientation, plane, unit cell | World-space text and shapes; source atom picking and PAE pair selection are provided |
| Computed properties | geometric interactions, Shrake–Rupley ASA, explicit cross-link bounds; optional overlap clashes | Geometry-based estimates differ from Mol* feature typing and validation services; annotation inputs are preferred for external evaluations |
| Structural extensions | membrane, assembly symmetry, confal pyramids, NtC tube, clashes, tunnels | Consume local results; no ANVIL optimization, symmetry server, DNATCO assignment, or tunnel search is inferred |
| Alpha orbitals | real spherical Gaussian basis L=0..4; orbital and occupied density surfaces | CPU collocation; same reference polynomials and order conventions; finite grid/cutoff/isovalue sampling |
| Meshes/segmentation | indexed JSON/NPZ and mesh CIF/BCIF; integer label grids | Local geometry only; no segmentation editing UI or streaming server |
| Kinemage/G3D | primitive list drawing and local compressed genome coordinates | No browser controls; Kinemage ribbons use tube segments |
| MVS | local structure/component/representation/color/opacity/transform/label trees and snapshot selection | Explicit subset of MVS protocol; annotations/complex selectors require normalized schemas; unsupported nodes raise errors |
| Quality/charge views | geometry, density, RCI, pLDDT, QMEAN, partial charges | Use supplied annotations; no quality/charge model is evaluated automatically |
| Automatic presets | auto, auto-lod, mesoscale | Static atom-count thresholds select local presets; no Mol* size categorization or camera-dependent LOD |
| Pairwise metrics | PAE/dense/sparse matrix panel, missing-data color, residue picking, PNG export | Separate Qt panel; one matrix export at a time |
| Themes | 39 built-in color names, 6 size names, extension annotation palettes; whole-object chain order grouped by mmCIF entity | Independent interpolation/Lab adjustment; absent mmCIF entity/order metadata cannot be recovered from PDB; themes require compatible data |
| Materials | Mol* GGX/Schlick/Smith, diffuse/metal ambient terms, camera-relative lights; GPU noise-based bumps and flat shading | Dedicated ray evaluates at vertices without fragment derivatives or bump; standard ray retains native lighting. Mixed unmanaged ray scenes retain native global illumination |
| Screen effects | Mol* 32-sample SSAO with bilateral blur; outline, shadow, cel, xray, unlit, bloom, DOF, illumination, antialias | Layer-local effects; no transparent or multiscale SSAO; direct-volume layers skip SSAO. Other effects and AA remain approximations; illumination is not a path tracer |
| Background | gradient, radial gradient, local image, six local skybox faces | Prepared bitmap approximation; skybox updates through the maintenance timer; dedicated ray export composites backgrounds |
| Clipping | six half-spaces | Shared uncapped clipped triangles and masked volume voxels; not a full clipping-object editor |

## Validation

`tests/reference` executes the pinned Mol* TypeScript in a local browser, exports
its geometry and camera, and compares actual PyMOL GPU/ray output. See the
[comparison report](fidelity.md) for conditions, metrics and portable commands.

`tests/test_reference_geometry.py` checks dimensions, tensor equations, normals,
signed fields, actual colors, PCA centers, angle-sector area, and target transforms.
`tests/test_geometry.py` covers every standard visual branch and every drawing
family. `tests/test_data_themes.py` checks theme outputs, orbital/ASA/measurement
formulas, map axes/units, BinaryCIF masks, and local MVS. `tests/test_lifecycle.py`
uses real PyMOL to verify overlaps, failures, manual deletion, sessions, and states.
`tests/check_real.py --gui` records GPU and ray images and fails on empty output or
OpenGL errors. `--visuals` tests the standard subvisuals separately. Regenerable
reports, temporary images, and benchmark results are deliberately excluded from Git.
Published gallery PNGs are tracked and regenerated locally, without CI.

The 500-residue, 100-state benchmark uses a repeated protein backbone with synthetic
displacements and a 1280×720 viewport. It reports preparation time, both CPU payload
sizes, bounded GPU cache use, rotation/state-switch FPS, actual movie callback state
rate, and peak RSS. Re-run on the target machine; frame rates are not universal.
