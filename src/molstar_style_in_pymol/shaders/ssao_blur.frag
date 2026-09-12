#version 120
// Mol* depth-aware Gaussian blur (MIT; see NOTICE).
varying vec2 uv;
uniform sampler2D depth;
uniform sampler2D ao;
uniform mat4 inverseProjection;
uniform vec2 pixel;
uniform vec2 direction;
uniform float depthBias;
uniform int kernelSize;
vec3 position(vec2 at,float z) {
    vec4 p=inverseProjection*vec4(at*2.0-1.0,z*2.0-1.0,1.0);
    return p.xyz/p.w;
}
void main() {
    float z=texture2D(depth,uv).r;
    if(z>=1.0 || z<=0.0) { gl_FragColor=vec4(1.0); return; }
    vec3 p=position(uv,z);
    float pixelSize=distance(p,position(uv+vec2(pixel.x,0.0),z));
    float sum=0.0,weight=0.0,sigma=float(kernelSize)/3.0;
    for(int i=-12;i<=12;i++) {
        float x=float(i);
        if(abs(x)>float(kernelSize/2) || (abs(x)>1.0 && abs(x)*pixelSize>0.8)) continue;
        vec2 at=uv+x*direction*pixel;
        if(at.x<0.0 || at.y<0.0 || at.x>1.0 || at.y>1.0) continue;
        float sd=texture2D(depth,at).r;
        if(sd>=1.0 || sd<=0.0 || abs(position(at,sd).z-p.z)>=depthBias) continue;
        float w=exp(-x*x/(2.0*sigma*sigma));
        sum+=texture2D(ao,at).r*w; weight+=w;
    }
    gl_FragColor=vec4(vec3(weight>0.0?sum/weight:1.0),1.0);
}
