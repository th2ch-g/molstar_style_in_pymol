#version 120
#include lighting
varying vec3 normal;
varying vec3 position;
varying vec3 world;
varying vec4 color;
uniform vec3 material;
uniform int unlit;
uniform int cel;
uniform float celSteps;
uniform int xray;
uniform int flatShaded;
uniform int clipCount;
uniform vec4 clipPlanes[6];
uniform int lightCount;
uniform vec3 lightDirection[8];
uniform vec3 lightColor[8];
uniform vec3 ambientColor;
uniform float exposure;
uniform float bumpFrequency;
uniform float bumpAmplitude;
void main() {
    for (int i=0; i<6; i++) {
        if (i < clipCount && dot(vec4(world, 1.0), clipPlanes[i]) < 0.0) discard;
    }
    vec3 n = normalize(normal);
    if (flatShaded == 1) n = normalize(cross(dFdx(position), dFdy(position)));
    else if (!gl_FrontFacing) n = -n;
    if (material.z > 0.0 && bumpFrequency > 0.0 && bumpAmplitude > 0.0)
        n = perturbNormal(-position, n, fbm(world*bumpFrequency), bumpAmplitude*material.z/bumpFrequency);
    vec3 v = normalize(-position);
    GeometricContext geometry;
    geometry.position = position;
    geometry.normal = n;
    geometry.viewDir = v;
    PhysicalMaterial physical;
    physical.diffuseColor = color.rgb * (1.0 - material.x);
    vec3 dxy = max(abs(dFdx(n)), abs(dFdy(n)));
    float geometryRoughness = max(max(dxy.x, dxy.y), dxy.z);
    physical.roughness = min(max(material.y, 0.0525) + geometryRoughness, 1.0);
    physical.specularColor = mix(vec3(0.04), color.rgb, material.x);
    physical.specularF90 = 1.0;
    ReflectedLight reflected = ReflectedLight(vec3(0.0), vec3(0.0), vec3(0.0), vec3(0.0));
    IncidentLight light;
    vec3 outgoing = vec3(0.0);
    for (int i=0; i<8; i++) {
        if (i >= lightCount) break;
        light.direction = lightDirection[i];
        light.color = lightColor[i] * PI;
        if (cel == 1) {
            float nl = max(dot(n, light.direction), 0.0);
            vec3 spec = nl * BRDF_GGX(light.direction, v, n, physical.specularColor, 1.0, max(material.y, 0.05));
            float intensity = nl * RECIPROCAL_PI * (1.0-material.x) + dot(spec, vec3(0.2126,0.7152,0.0722));
            outgoing += color.rgb * light.color * ceil(intensity * celSteps) / celSteps;
        } else {
            RE_Direct_Physical(light, geometry, physical, reflected);
        }
    }
    if (cel == 1) outgoing += physical.diffuseColor * ambientColor;
    else {
        RE_IndirectDiffuse_Physical(ambientColor * PI, geometry, physical, reflected);
        RE_IndirectSpecular_Physical(ambientColor * material.x, ambientColor * material.x, vec3(0.0), geometry, physical, reflected);
        outgoing = reflected.directDiffuse + reflected.indirectDiffuse + reflected.directSpecular + reflected.indirectSpecular;
    }
    outgoing = clamp(outgoing, 0.01, 0.99);
    if (unlit == 1) outgoing = color.rgb;
    if (xray == 1) outgoing = mix(color.rgb*0.15, color.rgb, pow(1.0-abs(dot(n,v)), 1.5));
    gl_FragColor = vec4(outgoing * exposure, color.a);
}
