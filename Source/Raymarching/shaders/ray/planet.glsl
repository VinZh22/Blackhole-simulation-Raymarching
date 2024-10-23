#version 330 core

//
////// Ray marching shader
// Implements a randomely generated planet

// Inputs coming from: 2d uv Coord
in vec2 uvCoord;

// Output of the fragment shader - output color
out vec4 fragColor;

// Uniform values sent by the CGP based c++ code
uniform float iTime;     // time since simulation start

// Uniform values sent by the camera model
uniform float aspect_ratio;
uniform float fov;
uniform float DIST_TO_CANVAS;
uniform vec3 ro;
uniform mat3 orientation;

const int MAX_MARCHING_STEPS = 256;
const float MIN_DIST = 0.0;
const float MAX_DIST = 1000.0;
const float PRECISION = 0.001;
const float dX = 0.001;

const float PI = 3.14159265359;

const int SEED = 0;

const vec3 SPHERE_CENTER = vec3(0,0, 0);
const float SPHERE_R = 1;
const float GROUND_HEIGHT = -15.0;
vec3 LIGHT_POSITION = vec3(45.0,15.0,45.0);
vec3 LIGHT_COLOR = vec3(1.0,1.0,1.0);

mat3 rot_y(float x) {
    return mat3(cos(x),0,sin(x),
                0,1,0,
                -sin(x),0, cos(x));
}

float rand(vec2 co) {
    return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

float rand_vals(vec2 x) {
  return rand(x + SEED);
}

// Need for better value noise function with splines
float value_noise_2D(vec2 x, int N, float scaling) { // takes an input in [0,1] x [0,1], returns value noise
  x = N*x;
  vec2 x_int = floor(x);
  vec2 x_fract = fract(x);
  float i = x_int.x;
  float j = x_int.y;
  
  return mix( mix(rand_vals(vec2(i,j)),
                  rand_vals(vec2(mod(i+1,N), j)),
                  x_fract.x),
              mix(rand_vals(vec2(i,mod(j+1,N))),
                  rand_vals(vec2(mod(i+1,N), mod(j+1,N))),
                  x_fract.x),
              x_fract.y) * scaling;
}

float displace(vec3 p) {
  vec3 direction = normalize(p);
  vec2 phi_theta = vec2(acos(direction.y), PI/2 + atan(direction.x/direction.z));
  if (direction.z < 0) phi_theta.y += PI;

  phi_theta /= vec2(PI, 2*PI);

  float displacement = value_noise_2D(phi_theta, 40, 1);
  displacement += value_noise_2D(phi_theta, 80, 0.5);
  displacement += value_noise_2D(phi_theta, 160, 0.25);
  displacement += value_noise_2D(phi_theta, 320, 0.125);
  displacement += value_noise_2D(phi_theta, 640, 0.0675);

  return displacement*sin(PI*phi_theta.x)*SPHERE_R*0.025;
}

float sdPlanet(vec3 p, vec3 center, float r ) {
  p = rot_y(0.1*iTime)*(p-center);

  float d1 = length(p) - r;
  float d2 = displace(p);

  return d1 + d2;
}

float sdfScene(vec3 p) {
  return sdPlanet(p, SPHERE_CENTER, SPHERE_R);
}

vec3 gradient(vec3 p) {
    vec3 dx = vec3(dX,0,0);
    vec3 dy = vec3(0,dX,0);
    vec3 dz = vec3(0,0,dX);
    vec3 grad;
    grad.x = (sdfScene(p+dx)-sdfScene(p))/dX;
    grad.y = (sdfScene(p+dy)-sdfScene(p))/dX;
    grad.z = (sdfScene(p+dz)-sdfScene(p))/dX;
    return grad;
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
    float d = min(  sdPlanet(p, SPHERE_CENTER, SPHERE_R),
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
  col = noise3(rd);

  if (d > MAX_DIST) {
    // The ray didn't hit anything: We display some background
  } else {
    // The ray hit the sphere or the ground, we need to compute the impact point
    vec3 p = ro + rd * d; // Hit point on some object
    if (sdPlanet(p ,SPHERE_CENTER, SPHERE_R) < 2*PRECISION) {
      // The ray hit the planet
      // Compute directions
      vec3 light_direction = normalize(LIGHT_POSITION-p);
      vec3 normal = normalize(gradient(p));
      vec3 reflect_direction = reflection(-rd, normal);

      // Compute phong reflexion and diffusion
      float diffus = clamp(dot(normal, light_direction),0,1);
      float reflected = pow(clamp(dot(reflect_direction,
                                      light_direction),
                                  0, 20),
                            24);
      
      // Compute shadows
      vec3 new_ro = p + normal * 2 * PRECISION;
      float shadow_d = rayMarch(new_ro, light_direction, MIN_DIST, MAX_DIST); // cast shadow ray to the light source
      if (shadow_d < length(LIGHT_POSITION - new_ro)) diffus *= 0.; // if the shadow ray hits the sphere, set
      
      col = vec3(0.5) * LIGHT_COLOR * clamp((diffus+reflected), 0.3, 1);

    } else {
      // The ray hit the ground
      // Compute directions
      vec3 light_direction = normalize(LIGHT_POSITION-p);
      vec3 normal = normalGround();
      vec3 reflect_direction = reflection(-rd, normal);

      // Compute phong diffusion
      float diffus = clamp(dot(normal, light_direction),0,1);
      
      // Compute shadows
      vec3 new_ro = p + normal * 2 * PRECISION;
      float shadow_d = rayMarch(new_ro, light_direction, MIN_DIST, MAX_DIST); // cast shadow ray to the light source
      if (shadow_d < length(LIGHT_POSITION - new_ro)) diffus *= 0.; // if the shadow ray hits the sphere, set
      
      col = LIGHT_COLOR * clamp((diffus), 0.3, 1);
    }
  }

  // Output to screen
  fragColor = vec4(col, 1.0);
}
