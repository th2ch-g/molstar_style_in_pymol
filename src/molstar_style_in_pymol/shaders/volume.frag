#version 120
varying vec2 uv;
uniform sampler3D field;
uniform sampler1D transfer;
uniform sampler2D sceneDepth;
uniform mat4 inverseMvp;
uniform mat4 mvp;
uniform mat4 worldToGrid;
uniform vec3 dimensions;
uniform vec2 domain;
uniform float stepSize;
uniform float opacity;
void main() {
    vec4 h0=inverseMvp*vec4(uv*2.0-1.0,-1.0,1.0);
    vec4 h1=inverseMvp*vec4(uv*2.0-1.0,1.0,1.0);
    vec3 a=h0.xyz/h0.w, b=h1.xyz/h1.w;
    vec3 direction=normalize(b-a);
    vec3 o=(worldToGrid*vec4(a,1.0)).xyz;
    vec3 d=(worldToGrid*vec4(direction,0.0)).xyz;
    vec3 safe=sign(d+vec3(1e-12))*max(abs(d),vec3(1e-9));
    vec3 t0=(vec3(0.0)-o)/safe, t1=(dimensions-1.0-o)/safe;
    vec3 near=min(t0,t1), far=max(t0,t1);
    float start=max(0.0,max(near.x,max(near.y,near.z)));
    float end=min(far.x,min(far.y,far.z));
    if (start>=end) discard;
    float oldZ=texture2D(sceneDepth,uv).r;
    vec4 occluder=inverseMvp*vec4(uv*2.0-1.0,oldZ*2.0-1.0,1.0);
    end=min(end,length(occluder.xyz/occluder.w-a));
    float delta=max(stepSize,(end-start)/1023.0);
    vec4 sum=vec4(0.0);
    float first=start;
    for (int i=0;i<1024;i++) {
        float t=start+float(i)*delta;
        if(t>end || sum.a>0.985) break;
        vec3 idx=o+t*d;
        float value=texture3D(field,(idx+0.5)/dimensions).r;
        vec4 rgba=texture1D(transfer,clamp((value-domain.x)/(domain.y-domain.x),0.0,1.0));
        rgba.a=(1.0-pow(1.0-rgba.a,delta/stepSize))*opacity;
        if(sum.a<0.01 && rgba.a>0.01) first=t;
        sum.rgb+=(1.0-sum.a)*rgba.a*rgba.rgb;
        sum.a+=(1.0-sum.a)*rgba.a;
    }
    if(sum.a<0.005) discard;
    vec4 clip=mvp*vec4(a+first*direction,1.0);
    gl_FragDepth=(clip.z/clip.w+1.0)*0.5;
    gl_FragColor=vec4(sum.rgb/max(sum.a,1e-6),sum.a);
}
