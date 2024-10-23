#version 330 core

//
////// Ray marching shader
// Implements star generation with perlin noise and relativistic bent light rays

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

const int MAX_MARCHING_STEPS = 255;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float MAX_DS = 0.1;
const float PRECISION = 0.001;

const float PI = 3.14159265359;

float BH_RADIUS = .5; // = to 2GM/c² the blackhole Schwarzschild radius


float rand(vec2 co) {
    return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

float sdSphere(vec3 p, vec3 center, float r ) {
  return length(p-center) - r;
}

float sdScene(vec3 p) {
  return sdSphere(p, vec3(0,0,0), BH_RADIUS);
}

float n (vec3 p) {
  return 1/(1-(BH_RADIUS-PRECISION)/length(p));
}

vec3 grad_n (vec3 p) {
  float n = n(p);
  float r = length(p);
  return -n*n*BH_RADIUS*p/r/r/r;
}

vec3[2] rayMarch(	vec3 ro, 		// ray origin
                vec3 rd, 		// ray direction
                float start,
                float end) {
  vec3 velocity = rd;
  vec3 p = ro;
  float ds = MAX_DS;

  for (int i = 0; i < MAX_MARCHING_STEPS; i++) { // dx = v dt, dv =  F dt, we want dx to be bounded by the sdf distance, we use variable dt
    float d = sdScene(p);
    if (d < PRECISION || length(p-ro) > end) break;

    float refractive_idx = n(p);
    ds = d/length(velocity/refractive_idx);
    ds = min(MAX_DS, ds);

    p += velocity*ds/refractive_idx;
    velocity += grad_n(p) * ds;
  }

  vec3[2] pos_dir;
  pos_dir[0] = p;
  pos_dir[1] = normalize(velocity);

  return pos_dir;
}

vec3 stars(
  vec3 rd, //ray direction
  int grid_num, //number of grid points for the stars
  float inverse_width, //The higher this float the smaller the stars
  float scarcity_factor, //The higher this float the more likely a star's intesity is to be reduced
  int inverse_band_width, //The higher this float narrower the generated star band will be
  float global_scaling, //A global scaling factor for all the stars intensity
  float seed = 1
) {
  float rand_angle_1 = 0.3*rand(vec2(inverse_width,scarcity_factor)+seed)-0.15;
  float c1 = cos(rand_angle_1);
  float s1 = sin(rand_angle_1);
  mat3 rot1 = mat3( 1,0,0,
                    0,c1,-s1,
                    0,s1,c1);

  float rand_angle_2 = 0.3*rand(vec2(inverse_width,scarcity_factor)+seed)-0.15;
  float c2 = cos(rand_angle_2);
  float s2 = sin(rand_angle_2);
  mat3 rot2 = mat3( c2,-s2,0,
                    s2,c2,0,
                    0,0,1);

  vec3 rd_rot = rot2 * rot1 * rd;
  vec2 phi_theta = vec2(acos(rd_rot.y), atan(rd_rot.x/rd_rot.z));
  vec2 phi_theta_mod = mod(grid_num*phi_theta, PI);
  // generate stars at all grid points
  vec3 stars =  vec3( exp(0.2-inverse_width*length(phi_theta_mod-vec2(PI/2)))
                      * pow(sin(phi_theta.x), inverse_band_width));
  // randomely remove some stars
  stars *= global_scaling*exp(-scarcity_factor*rand(floor(grid_num*phi_theta/PI)+seed));
  // add random redshift for some color
  float random_factor = 2*rand(floor(1+(grid_num)*phi_theta/PI)+seed) - 1;
  stars *= vec3(1+0.3*random_factor,0.9,1-0.3*random_factor);

  return stars;
}

vec3 create_stars (vec3 rd) {
  vec3 col = stars(rd,10,9.0,8.0,4,2.0);
  col += stars(rd,20,9.0,10,6,1.0);
  col += stars(rd,20,9.0,10,6,1.0,2);
  col += stars(rd,40,12.0,12,6,1.0);
  col += stars(rd,40,12.0,12,6,1.0,2);
  col += stars(rd,40,12.0,12,6,1.0,3);
  col += stars(rd,40,12.0,12,6,1.0,4);
  col += stars(rd,80,12.0,12,8,1.0);
  col += stars(rd,80,12.0,12,8,1.0,2);
  col += stars(rd,80,12.0,12,8,1.0,3);
  col += stars(rd,80,12.0,12,8,1.0,4);
  col += stars(rd,120,12.0,14,8,1.0);
  col += stars(rd,120,12.0,14,8,1.0,2);

  return col;
}

void main() {
  vec2 uv = uvCoord;
  
  // We transform the coordinates so they fit the specified aspect ratio and fov
  vec3 uv3d = vec3(uv.x, uv.y, -DIST_TO_CANVAS);
  vec3 col = vec3(0); // Output color vector
  // vec3 rd = normalize(vec3(0.1*uv3d.x, 0.1*uv3d.y/aspect_ratio, uv3d.z)); // ray direction
  vec3 rd = normalize(vec3(fov*uv3d.x, fov*uv3d.y/aspect_ratio, uv3d.z)); // ray direction
  rd = normalize(orientation * rd);

  vec3[2] pos_dir = rayMarch(ro, rd, MIN_DIST, MAX_DIST); // position and direction of ray at end of raymarch
  float d = sdScene(pos_dir[0]); // distance to sdf

  if (d > PRECISION) {
    // The ray didn't hit anything: We display stars at infinity.

    // Galaxy-like background
    vec3 galaxy_bg = vec3(0.1*pow(cos(pos_dir[1].y),32));

    // Randomely generated stars through a perlin noise method.
    vec3 rand_stars = create_stars(pos_dir[1]);

    col = rand_stars+galaxy_bg;

  } else { 
    // The ray hit the black hole, we display pure black
    col = vec3(0);
  }

  // Output to screen
  fragColor = vec4(col, 1.0);
}
