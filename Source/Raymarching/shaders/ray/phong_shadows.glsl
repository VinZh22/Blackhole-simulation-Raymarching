#version 330 core

//
////// Ray marching shader
// Implements shading, phong reflection and shadows

// Inputs coming from: 2d uv Coord
in vec2 uvCoord;

// Output of the fragment shader - output color
out vec4 fragColor;

// Uniform values sent by the CGP based c++ code
uniform float time;     // time since simulation start

// Uniform values sent by the camera model
uniform float aspect_ratio;
uniform float fov;
uniform float DIST_TO_CANVAS;
uniform vec3 ro;
uniform mat3 orientation;

const int MAX_MARCHING_STEPS = 255;
const float MIN_DIST = 0.0;
const float MAX_DIST = 1000.0;
const float PRECISION = 0.01;

const float PI = 3.14159265359;

vec3 SPHERE_CENTER = vec3(0,0, 0);
float SPHERE_R = 1.0;
float GROUND_HEIGHT = -3.0;
vec3 LIGHT_POSITION = vec3(0, 4, 2);
vec3 LIGHT_COLOR = vec3(1.0,1.0,1.0);


float sdSphere(vec3 p, vec3 center, float r ) {
  return length(p-center) - r;
}

vec3 normalSphere(vec3 p) {
  return normalize(p-SPHERE_CENTER);
}

float sdGround(vec3 p, float height) {
  return length(p.y-height);
}

vec3 normalGround() {
  return vec3(0,1,0);
}

float rayMarch(	vec3 ro, 		// ray origin
                vec3 rd, 		// ray direction
                float start,
                float end) {
  float depth = start;

  for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
    vec3 p = ro + depth * rd;
    float d = min(  sdSphere(p, SPHERE_CENTER, SPHERE_R),
                    sdGround(p, GROUND_HEIGHT));
    depth += d;
    if (d < PRECISION || depth > end) break;
  }

  return depth;
}

vec3 reflection(vec3 direction, vec3 normal) {
  vec3 direction_n = normal*dot(direction, normal);
  return 2*direction_n - direction;
}

void main() {
  vec2 uv = uvCoord;
  
  // We transform the coordinates so they fit the specified aspect ratio and fov
  vec3 uv3d = vec3(uv.x, uv.y, -DIST_TO_CANVAS);
  vec3 col = vec3(0); // Output color vector
  vec3 rd = normalize(vec3(fov*uv3d.x, fov*uv3d.y/aspect_ratio, uv3d.z)); // ray direction
  rd = normalize(orientation * rd);

  float d = rayMarch(ro, rd, MIN_DIST, MAX_DIST); // distance to sphere

  // By default we display some blue sky background
  col = vec3(0);

  if (d > MAX_DIST) {
    // The ray didn't hit anything: We display some background
  } else {
    // The ray hit the sphere or the ground, we need to compute the impact point
    vec3 p = ro + rd * d; // Hit point on some object
    if (sdSphere(p ,SPHERE_CENTER, SPHERE_R) < 2*PRECISION) {
      // The ray hit the sphere
      // Compute directions
      vec3 light_direction = normalize(LIGHT_POSITION-p);
      vec3 normal = normalSphere(p);
      vec3 reflect_direction = reflection(-rd, normal);

      // Compute phong reflexion and diffusion
      float dist = length(p-LIGHT_POSITION);
      float diffus = clamp(dot(normal, light_direction),0,1);
      float reflected = pow(clamp(dot(reflect_direction,
                                          light_direction),
                                      0, 20),
                                24);
      
      // Compute shadows
      vec3 new_ro = p + normal * 2 * PRECISION;
      float shadow_d = rayMarch(new_ro, light_direction, MIN_DIST, MAX_DIST); // cast shadow ray to the light source
      if (shadow_d < length(LIGHT_POSITION - new_ro)) diffus *= 0.; // if the shadow ray hits the sphere, set

      // Distance based attenuation
      float attenuation = clamp(pow(dist/100,-2), 0, 1);
      
      col = vec3(1.0,0.2,0.2) * LIGHT_COLOR * clamp((diffus+reflected)*attenuation, 0.3, 1);

    } else {
      // The ray hit the ground
      // Compute directions
      vec3 light_direction = normalize(LIGHT_POSITION-p);
      vec3 normal = normalGround();
      vec3 reflect_direction = reflection(-rd, normal);

      // Compute phong diffusion
      float dist = length(p-LIGHT_POSITION);
      float diffus = clamp(dot(normal, light_direction),0,1);
      
      // Compute shadows
      vec3 new_ro = p + normal * 2 * PRECISION;
      float shadow_d = rayMarch(new_ro, light_direction, MIN_DIST, MAX_DIST); // cast shadow ray to the light source
      if (shadow_d < length(LIGHT_POSITION - new_ro)) diffus *= 0.; // if the shadow ray hits the sphere, set

      // Distance based attenuation
      float attenuation = clamp(pow(dist/100,-2), 0, 1);
      
      col = LIGHT_COLOR * clamp((diffus)*attenuation, 0.3, 1);
    }
  }

  // Output to screen
  fragColor = vec4(col, 1.0);
}
