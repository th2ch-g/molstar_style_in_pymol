#version 120
varying vec3 normal;
varying vec3 position;
varying vec3 world;
varying vec4 color;
uniform vec3 material;
uniform int unlit;
uniform int cel;
uniform int xray;
uniform int flatShaded;
uniform int clipCount;
uniform vec4 clipPlanes[6];
void main() {
    for (int i=0; i<6; i++) {
        if (i < clipCount && dot(vec4(world, 1.0), clipPlanes[i]) < 0.0) discard;
    }
    vec3 n = normalize(normal);
    if (!gl_FrontFacing) n = -n;
    if (material.z > 0.0) n = normalize(n + material.z * 0.1 * sin(world * 14.0));
    vec3 v = normalize(-position);
    vec3 l = normalize(vec3(-0.35, 0.6, 1.0));
    float diffuse = max(0.0, dot(n, l));
    if (cel == 1) diffuse = floor(diffuse * 4.0) / 3.0;
    float specular = pow(max(0.0, dot(n, normalize(v+l))), mix(140.0, 6.0, material.y));
    vec3 base = color.rgb;
    vec3 outColor = base * (0.28 + 0.72 * diffuse) * (1.0 - material.x*0.35);
    outColor += mix(vec3(1.0), base, material.x) * specular * (1.0-material.y*0.65) * 0.65;
    if (unlit == 1) outColor = base;
    if (xray == 1) outColor = mix(base*0.15, base, pow(1.0-abs(dot(n,v)), 1.5));
    gl_FragColor = vec4(clamp(outColor, 0.0, 1.0), color.a);
}
