#version 120
varying vec2 uv;
uniform sampler2D image;
uniform sampler2D depth;
uniform sampler2D ambientOcclusion;
uniform vec2 pixel;
uniform float occlusion;
uniform float outline;
uniform float shadow;
uniform float bloom;
uniform float dof;
uniform float focus;
uniform float antialias;
uniform float illumination;
void main() {
    vec4 c = texture2D(image, uv);
    float z = texture2D(depth, uv).r;
    if (c.a < 0.001) discard;
    c.rgb /= c.a;
    float ao = 0.0;
    float edge = 0.0;
    vec4 blur = vec4(0.0);
    vec3 glow = vec3(0.0);
    for (int i=0; i<12; i++) {
        float a = float(i) * 6.2831853 / 12.0;
        vec2 delta = vec2(cos(a), sin(a)) * pixel;
        float d = texture2D(depth, uv + delta*3.0).r;
        float difference=z-d;
        ao += smoothstep(0.00008,0.002,difference)*(1.0-smoothstep(0.01,0.025,difference));
        edge += abs(d-z) > 0.003 ? 1.0 : 0.0;
        vec4 sampleColor = texture2D(image, uv + delta*max(1.0, abs(z-focus)*dof*200.0));
        blur += sampleColor;
        vec3 light = texture2D(image, uv + delta*5.0).rgb;
        glow += light * max(0.0, max(light.r, max(light.g, light.b))-0.65);
    }
    c.rgb *= max(0.01, 1.0-occlusion*(1.0-texture2D(ambientOcclusion,uv).r));
    c.rgb *= 1.0-min(0.8,illumination*ao/32.0);
    c.rgb *= 1.0 - outline * min(1.0,edge/4.0);
    float shade = texture2D(depth, uv + pixel*vec2(-8.0,12.0)).r;
    c.rgb *= shade<z-0.0005 && shade>z-0.08 ? 1.0-shadow*0.3 : 1.0;
    c.rgb += glow*bloom/12.0;
    if (dof > 0.0) c = mix(c, blur/12.0, min(0.85,abs(z-focus)*dof*10.0));
    if (antialias>0.0 && edge>0.0) c=mix(c,blur/12.0,0.2);
    gl_FragColor=c;
    gl_FragDepth=z;
}
