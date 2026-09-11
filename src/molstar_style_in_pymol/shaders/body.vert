#version 120
varying vec3 normal;
varying vec3 position;
varying vec3 world;
varying vec4 color;
void main() {
    world = gl_Vertex.xyz;
    position = (gl_ModelViewMatrix * gl_Vertex).xyz;
    normal = normalize(gl_NormalMatrix * gl_Normal);
    color = gl_Color;
    gl_Position = ftransform();
}
