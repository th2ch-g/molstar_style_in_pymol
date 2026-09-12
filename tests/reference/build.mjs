// Bundle a local, pinned Mol* checkout for development comparisons only.
import { createRequire } from 'node:module';
import { resolve } from 'node:path';
import { mkdir, writeFile } from 'node:fs/promises';
import { execFileSync } from 'node:child_process';

const root = resolve(process.argv[2] || '../molstar');
const output = resolve(process.argv[3] || '.cache/reference');
const require = createRequire(`${output}/package.json`);
const revision = execFileSync('git', ['rev-parse', 'HEAD'], { cwd: root, encoding: 'utf8' }).trim();
if (revision !== '5b1b54ed03b03936041f514b33b8bb129b774d37') throw new Error('Use the documented Mol* reference revision');
const { build } = require('esbuild');
await mkdir(output, { recursive: true });
await build({
    entryPoints: ['tests/reference/viewer.ts'],
    outfile: `${output}/viewer.js`,
    bundle: true,
    platform: 'browser',
    nodePaths: [`${output}/node_modules`],
    alias: { 'molstar-reference': `${root}/src` },
    define: { 'process.env.NODE_ENV': '"production"', '__REFERENCE_REVISION__': JSON.stringify(revision) },
    logLevel: 'warning',
});
await writeFile(`${output}/index.html`, `<!doctype html><html><head><meta charset="utf-8">
<style>html,body,#viewer{margin:0;width:100%;height:100%;overflow:hidden}canvas{display:block;width:100%;height:100%}</style>
</head><body><div id="viewer"><canvas id="canvas"></canvas></div><script src="viewer.js"></script></body></html>`);
