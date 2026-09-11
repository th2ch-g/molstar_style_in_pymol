#version 120
varying vec2 uv;
uniform sampler2D image;
void main() { gl_FragColor=texture2D(image,uv); gl_FragDepth=1.0; }
