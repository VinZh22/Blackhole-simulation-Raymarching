#version 330 core

in vec2 fragCoord;
//uniform vec2 iMouse;
uniform vec2 iResolution;
uniform float aspect_ratio;
uniform float iTime;

out vec4 fragColor;

//vec2 fragCoord = gl_FragCoord.xy;

#define PI 3.14159265359
#define TWOPI 6.28318530718
//drag the window LR to control roughness

//--graphics setting (lower = better fps)---------------------------------------------------------------------
const int MAX_MARCHING_STEPS = 255;
const float MIN_DIST = 0.0;
const float MAX_DIST = 50.0;
const float PRECISION = 0.001;

const float radius = 0.5;

float sdSphere(vec3 p, float r )
{
  vec3 offset = vec3(0, 0, -2);
  return length(p - offset) - r;
}

float rayMarch(vec3 ro, vec3 rd, float start, float end) {
    float depth = start;

    for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
        vec3 p = ro + depth * rd;
        float d = sdSphere(p, radius);
        depth += d;
        if (d < PRECISION || depth > end) break;
    }

    return depth;
}

vec3 calcNormal(vec3 p) {
    vec2 e = vec2(1.0, -1.0) * 0.0005; // epsilon
    float r = radius; // radius of sphere
    return normalize(
      e.xyy * sdSphere(p + e.xyy, r) +
      e.yyx * sdSphere(p + e.yyx, r) +
      e.yxy * sdSphere(p + e.yxy, r) +
      e.xxx * sdSphere(p + e.xxx, r));
}

void main()
{
    vec2 uv = fragCoord;//(fragCoord-.5*iResolution.xy)/iResolution.y;
    //uv = uv / vec2(1.0, aspect_ratio);
    vec3 backgroundColor = vec3(0., 1, 1);

    vec3 col = vec3(0);
    vec3 ro = vec3(0, 0, 1); // ray origin that represents camera position
    vec3 rd = normalize(vec3(uv, -1)); // ray direction

    float d = rayMarch(ro, rd, MIN_DIST, MAX_DIST); // distance to sphere
    
    if (d>0){
        col = normalize(vec3(0., d/50., 0.0));
    }

    else if (d > MAX_DIST) {
        col = backgroundColor; // ray didn't hit anything
    } else {
        vec3 p = ro + rd * d; // point on sphere we discovered from ray marching
        vec3 normal = calcNormal(p);
        vec3 lightPosition = vec3(2, 2, 7);
        vec3 lightDirection = normalize(lightPosition - p);

        // Calculate diffuse reflection by taking the dot product of
        // the normal and the light direction.
        float dif = clamp(dot(normal, lightDirection), 0.3, 1.);

        // Multiply the diffuse reflection value by an orange color and add a bit
        // of the background color to the sphere to blend it more with the background.
        col = dif * vec3(1, 0., 0.) + backgroundColor * .2;
        col = vec3(0., 1.0, 0.0);
    }

    // Output to screen
    fragColor = vec4(col, 1.0);
}
