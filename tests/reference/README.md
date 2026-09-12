# Local Mol* reference comparison

This development harness executes the original Mol* TypeScript. It is not a runtime
dependency and is not a CI workflow. Start from the repository root after `pixi install`.
Use a sibling `../molstar` checkout at revision
`5b1b54ed03b03936041f514b33b8bb129b774d37`; `build.mjs` rejects another revision.

Place RCSB mmCIF inputs `8gng.cif`, `1crn.cif` and `1bna.cif` under `.cache/reference/`.
They are available from the corresponding `files.rcsb.org/download/<ID>.cif` URLs.
These downloaded inputs and all raw reports are ignored.

```sh
npm install --prefix .cache/reference --no-audit --no-fund esbuild rxjs tslib io-ts fp-ts mutative immutable
node tests/reference/build.mjs
export PLAYWRIGHT_BROWSERS_PATH=.cache/reference/browsers
pixi run uv run --no-project --with playwright playwright install chromium
pixi run uv run --no-project python -m http.server 8765 --bind 127.0.0.1 --directory .cache/reference
```

With that local server running, use another terminal with the same browser path:

```sh
export PLAYWRIGHT_BROWSERS_PATH=.cache/reference/browsers
pixi run uv run --no-project --with playwright python tests/reference/render.py --suite tests/reference/cases.json
pixi run uv run --no-project python tests/reference/compare.py --structure .cache/reference/8gng.cif --reference .cache/reference/molstar-8gng-color.json
pixi run uv run --no-project python tests/reference/pymol_render.py --structure .cache/reference/8gng.cif --reference .cache/reference/molstar-8gng-color.json --occlusion --name pymol-8gng-color
pixi run uv run --no-project python tests/reference/metrics.py .cache/reference/molstar-8gng-color.png .cache/reference/pymol-8gng-color-gpu.png
```

`render.py --browser` accepts an installed Chromium executable as an alternative
to Playwright's browser download. The browser runs headlessly with SwiftShader;
PyMOL uses an independent actual Qt/OpenGL process. Never reuse a user's session.
Repeat the last three commands for the other named cases, including the `-ray.png`
output. All examples use SSAO, so pass `--occlusion`. `--native-ss` checks native
PyMOL assignments instead of applying the exported reference assignments.

`compare.py` exits nonzero when curve/frame or surface-vertex differences exceed
0.0001; cap triangulation centers are excluded from the surface test. `metrics.py`
reports registered-image silhouette IoU and RGB MAE, without fitting any image
transform or color correction. Its white-background mask and two-pixel interior
erosion are explicit in the source. Published PNGs and curated results are in
`docs/gallery/` and `docs/fidelity.md`; keep numerical exports under `.cache/`.
