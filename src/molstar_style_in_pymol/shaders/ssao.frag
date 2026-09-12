#version 120
// Mol* view-space SSAO, adapted from ssao.frag.ts (MIT; see NOTICE).
varying vec2 uv;
uniform sampler2D depth;
uniform vec2 pixel;
uniform mat4 projection;
uniform mat4 inverseProjection;
uniform vec3 samples[32];
uniform float radius;
uniform float bias;
vec3 viewPosition(vec2 at, float z) {
    vec4 p = inverseProjection * vec4(at*2.0-1.0,z*2.0-1.0,1.0);
    return p.xyz/p.w;
}
float readDepth(vec2 at) { return texture2D(depth,at).r; }
vec3 positionAt(vec2 at) { return viewPosition(at,readDepth(at)); }
float noise(vec2 at) { return abs(fract(sin(mod(dot(at,vec2(12.9898,78.233)),3.141592653589793))*43758.5453)); }
float smoother(float x) { x=clamp(x,0.0,1.0); return x*x*x*(x*(x*6.0-15.0)+10.0); }
void main() {
    float z=readDepth(uv);
    if(z>=1.0) { gl_FragColor=vec4(1.0); return; }
    vec3 p=viewPosition(uv,z);
    vec2 dx=vec2(pixel.x,0.0),dy=vec2(0.0,pixel.y);
    vec3 l=p-positionAt(uv-dx),r=positionAt(uv+dx)-p;
    vec3 d=p-positionAt(uv-dy),u=positionAt(uv+dy)-p;
    vec2 he=abs(vec2(2.0*readDepth(uv-dx)-readDepth(uv-2.0*dx),2.0*readDepth(uv+dx)-readDepth(uv+2.0*dx))-z);
    vec2 ve=abs(vec2(2.0*readDepth(uv-dy)-readDepth(uv-2.0*dy),2.0*readDepth(uv+dy)-readDepth(uv+2.0*dy))-z);
    vec3 n=normalize(cross(he.x<he.y?l:r,ve.x<ve.y?d:u));
    vec3 random=normalize(vec3(vec2(noise(uv),noise(uv+vec2(3.141592653589793,2.71828)))*2.0-1.0,0.0));
    vec3 tangent=normalize(random-n*dot(random,n));
    mat3 frame=mat3(tangent,cross(n,tangent),n);
    float sum=0.0,count=32.0;
    for(int i=0;i<32;i++) {
        vec3 s=p+frame*samples[i]*radius;
        vec4 projected=projection*vec4(s,1.0);
        vec2 at=projected.xy/projected.w*0.5+0.5;
        if(at.x<0.0 || at.y<0.0 || at.x>1.0 || at.y>1.0) { count-=1.0; continue; }
        float sd=readDepth(at);
        if(sd<1.0) {
            float vz=viewPosition(at,sd).z;
            sum+=step(s.z+0.025,vz)*smoother(radius/max(abs(p.z-vz),0.000001));
        }
    }
    gl_FragColor=vec4(vec3(clamp(1.0-bias*sum/max(count,1.0),0.01,1.0)),1.0);
}
