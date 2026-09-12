// Execute the original Mol* implementation, never a Python test double.
import { PluginContext } from 'molstar-reference/mol-plugin/context';
import { DefaultPluginSpec } from 'molstar-reference/mol-plugin/spec';
import { PolymerTraceIterator } from 'molstar-reference/mol-repr/structure/visual/util/polymer/trace-iterator';
import { createCurveSegmentState, interpolateCurveSegment } from 'molstar-reference/mol-repr/structure/visual/util/polymer/curve-segment';
import { SecondaryStructureType, isNucleic } from 'molstar-reference/mol-model/structure/model/types';
import { MolScriptBuilder as MS } from 'molstar-reference/mol-script/language/builder';
import { PostprocessingParams } from 'molstar-reference/mol-canvas3d/passes/postprocessing';
import { ParamDefinition as PD } from 'molstar-reference/mol-util/param-definition';

const plugin = new PluginContext(DefaultPluginSpec());
const ready = (async () => {
    await plugin.init();
    await plugin.initViewerAsync(document.querySelector('canvas')!, document.querySelector('#viewer')!);
    plugin.canvas3d!.setProps({
        renderer: { backgroundColor: 0xffffff },
        camera: { mode: 'orthographic', helper: { axes: { name: 'off', params: {} } } },
        cameraFog: { name: 'off', params: {} },
        postprocessing: { occlusion: { name: 'off', params: {} }, antialiasing: { name: 'off', params: {} } },
    } as any);
})();

async function load(options: any = {}) {
    await ready;
    await plugin.clear();
    plugin.canvas3d!.setProps({ postprocessing: {
        occlusion: { name: 'off', params: {} }, antialiasing: { name: 'off', params: {} },
    } } as any);
    const text = await (await fetch(options.file || '8gng.cif')).text();
    const data = await plugin.builders.data.rawData({ data: text });
    const trajectory = await plugin.builders.structure.parseTrajectory(data, 'mmcif');
    const model = await plugin.builders.structure.createModel(trajectory);
    const structure = await plugin.builders.structure.createStructure(model, { name: 'model', params: {} });
    let component = await plugin.builders.structure.tryCreateComponentStatic(structure, 'polymer');
    if (options.chains) {
        component = await plugin.builders.structure.tryCreateComponentFromExpression(component!, MS.struct.generator.atomGroups({
            'chain-test': MS.core.set.has([MS.set(...options.chains), MS.struct.atomProperty.macromolecular.auth_asym_id()]),
        }), 'reference-chains');
    }
    const typeParams = { quality: 'custom', linearSegments: 8, radialSegments: 16,
        visuals: ['polymer-trace'], ...options.params };
    const representation = await plugin.builders.structure.representation.addRepresentation(component!, {
        type: options.type || 'cartoon', typeParams,
        color: options.color || 'uniform', colorParams: { value: 0x1b9e77 },
        size: 'uniform', sizeParams: { value: 1 },
    });
    if (options.canvas) {
        for (const [key, value] of Object.entries(options.canvas.postprocessing || {}) as any) {
            const definition = (PostprocessingParams as any)[key];
            if (definition?.type === 'mapped') value.params = { ...PD.getDefaultValues(definition.map(value.name).params), ...value.params };
        }
        plugin.canvas3d!.setProps(options.canvas);
    }
    plugin.canvas3d!.requestCameraReset({ durationMs: 0 });
    plugin.canvas3d!.commit(true);
    if (options.camera) plugin.canvas3d!.camera.setState(options.camera, 0);
    (window as any).representation = representation;
    (window as any).structure = component!.obj!.data;
    (window as any).referenceOptions = options;
    return { camera: plugin.canvas3d!.camera.getSnapshot(), params: plugin.canvas3d!.props };
}

function dump() {
    const structure = (window as any).structure;
    const repr = (window as any).representation.obj.data.repr;
    const tubular = !!(window as any).referenceOptions.params?.tubularHelices;
    const meshes = repr.renderObjects.filter((o: any) => o.type === 'mesh').map((o: any) => ({
        vertices: Array.from(o.values.aPosition.ref.value), normals: Array.from(o.values.aNormal.ref.value),
        faces: Array.from(o.values.elements.ref.value), groups: Array.from(o.values.aGroup.ref.value),
        drawCount: o.values.drawCount.ref.value, vertexCount: o.values.uVertexCount.ref.value,
    }));
    const segments: any[] = [];
    for (const unit of structure.units) {
        if (unit.kind !== 0) continue;
        const it = PolymerTraceIterator(unit, structure, { ignoreSecondaryStructure: false, useHelixOrientation: tubular });
        while (it.hasNext) {
            const v = it.move();
            const h = unit.model.atomicHierarchy;
            const ri = h.residueAtomSegments.index[v.center.element];
            const ci = h.chainAtomSegments.index[v.center.element];
            const state = createCurveSegmentState(8);
            const helix = SecondaryStructureType.is(v.secStrucType, SecondaryStructureType.Flag.Helix);
            const sheet = SecondaryStructureType.is(v.secStrucType, SecondaryStructureType.Flag.Beta);
            interpolateCurveSegment(state, v, helix && !tubular ? 0.9 : 0.5, isNucleic(v.moleculeType) ? 0.3 : 0.5);
            segments.push({ chain: h.chains.auth_asym_id.value(ci), resi: String(h.residues.auth_seq_id.value(ri)) + h.residues.pdbx_PDB_ins_code.value(ri),
                kind: isNucleic(v.moleculeType) ? 'nucleic' : 'protein', trace: h.atoms.label_atom_id.value(v.center.element),
                ss: helix ? 'H' : sheet ? 'S' : '',
                controls: Object.fromEntries(['p0','p1','p2','p3','p4','d12','d23','secStrucFirst','secStrucLast','first','last','initial','final'].map(k => [k, Array.isArray(v[k]) ? [...v[k]] : v[k]])),
                curve: Array.from(state.curvePoints), normals: Array.from(state.normalVectors), binormals: Array.from(state.binormalVectors),
            });
        }
    }
    const samples = (plugin.canvas3dContext as any)?.passes?.draw?.postprocessing?.ssao?.renderable?.values?.uSamples?.ref?.value;
    return { revision: __REFERENCE_REVISION__, meshes, segments, canvas: plugin.canvas3d!.props, ssaoSamples: samples ? Array.from(samples) : null, camera: plugin.canvas3d!.camera.getSnapshot(), view: Array.from(plugin.canvas3d!.camera.view), projection: Array.from(plugin.canvas3d!.camera.projection) };
}
(window as any).reference = { ready, load, dump, plugin };
