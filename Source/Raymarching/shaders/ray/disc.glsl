#version 330 core

//
////// Ray marching shader
// Implements an acretion disc with volumetric cloud, perlin noise and animations.

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
const float PRECISION = 0.001;
const float MAX_DX = 0.1;

const float PI = 3.14159265359;

vec3 BH_CENTER = vec3(0.);
float BH_RAD = 1.0;

const float DISC_W_IN = 1.5;
const float DISC_W_OUT = 3.0;
const float DISC_H = .0;
const float RING_SPEED = 1.0;

const mat3 ID3 = mat3(1,0,0,
                            0,1,0,
                            0,0,1);


###
###### Sdf and scene definition
###

float sdSphere(vec3 p, vec3 center, float r ) {
  return length(p-center) - r;
}

float sdfTorus( vec3 p, vec2 t, vec3 center, mat3 orientation)
{
  p = orientation*(p-center);
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

float sdfCappedCylinder( vec3 p, float h, float r )
{
  vec2 d = abs(vec2(length(p.xz),p.y)) - vec2(r,h);
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdfBH (vec3 p) {
  return sdSphere(p, BH_CENTER, BH_RAD);
}

float sdfDisc (vec3 p) {
  return max(
    - sdfCappedCylinder(p, 1., DISC_W_IN),
    sdfCappedCylinder(p, DISC_H, DISC_W_OUT)
  );
}

float sdfScene(vec3 p) { // The scene contains all solid objects
  return min(
    sdfBH(p),
    sdfDisc(p)
  );
}



###
###### Cloud generation
###

vec4 cloud_torus(vec3 p) { // a torus shaped cloud using an  sdf
  p = p*vec3(1,2,1);
  float d = sdfTorus(
    p,
    vec2(DISC_W_IN+0.5,0.5),
    vec3(0,0,0),
    ID3
    );
  return vec4(0.8,0.8,1,0.5)*max(-d,0);
}

vec4 cloud_sphere(vec3 p) {
  float d = sdfBH(p);
  return vec4(0.8,0.8,1,0.5)*(0.5/(1+10*d*d));
}

vec4 clouds(vec3 p) {
  return cloud_torus(p) + cloud_sphere(p);
}



###
###### Acretion disc generation
###

float rand(vec2 co) {
    return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

//Second try: trace a very simple disc and substract lines afterwards
float lines_neg(
  float r,
  float theta,
  int num,
  float inverse_width, //The higher this float the smaller the line
  float displacement,
  float scarcity_factor, //The higher this float the more likely a star's intesity is to be reduced
  float global_scaling, //A global scaling factor for all the lines' intensity
  float seed
) {
  float delta = DISC_W_OUT-DISC_W_IN;
  float mod_r = mod(num*r, delta); // mod(x,y) = x - y * floor(x/y)
  float floor_r = floor(r*num)/num;

  // First we compute the opacity as an inverse of the distance of r to a grid on the disc making regularly spaced rings
  // We add a small random displacement to make the lines look less rigid
  float opacity = exp(-inverse_width * length(mod_r - delta/2 + displacement*rand(vec2(floor(num*r/delta) + seed))));
  
  // We randomely scale down the opacity of the rings to give them variable intensity
  opacity *= global_scaling*exp(-scarcity_factor*rand(vec2(floor(num*r/delta) + seed + 0.1)));

  // We scale down the opacity based on a random angle
  // We make that angle spin at a varying rate repending on the distance to the black hole
  float angle = 2*PI*rand(vec2(floor_r + seed + 0.2)) * iTime * RING_SPEED /(floor_r - 1) /(floor_r - 1);
  opacity*= min(exp(-length(mod(theta-angle,2*PI) - PI)), 1);

  return opacity;
}

vec4 disc_color(vec3 p) {
  float r = length(p.xz);
  float theta = atan(p.x/p.z);
  if (p.z < 0) theta += PI;

  vec4 color = vec4(.8,.8,1.,1);

  color.w -= lines_neg(r, theta, 10, 5, 1, 5, 5, 1);
  color.w -= lines_neg(r, theta, 12, 7, 1, 5, 5, 1);
  color.w -= lines_neg(r, theta, 12, 7, 1, 5, 5, 2);
  color.w -= lines_neg(r, theta, 14, 9, 1, 5, 5, 1);
  color.w -= lines_neg(r, theta, 16, 10, 1, 5, 5, 1);

  //We make the colors fade out on the edge of the disc
  color.w *= (1-exp(3*(r-DISC_W_OUT))) * (1-exp(-10*(r-DISC_W_IN)));

  //We alter the colors to absorb light instead of emmiting light on the edges
  float end = -0.2;
  color.xyz *= 1-(1-end)*(r-DISC_W_IN)/(DISC_W_OUT-DISC_W_IN);
  return color;
}



###
###### Random stars generation
###

vec3 stars(
  vec3 rd, //ray direction
  int grid_num, //number of grid points for the stars
  float inverse_width, //The higher this float the smaller the stars
  float scarcity_factor, //The higher this float the more likely a star's intesity is to be reduced
  int inverse_band_width, //The higher this float narrower the generated star band will be
  float global_scaling, //A global scaling factor for all the stars intensity
  float seed
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
  vec3 col = stars(rd,10,9.0,8.0,4,2.0,1);
  col += stars(rd,20,9.0,10,6,1.0,1);
  col += stars(rd,20,9.0,10,6,1.0,2);
  col += stars(rd,40,12.0,12,6,1.0,1);
  col += stars(rd,40,12.0,12,6,1.0,2);
  col += stars(rd,40,12.0,12,6,1.0,3);
  col += stars(rd,40,12.0,12,6,1.0,4);
  col += stars(rd,80,12.0,12,8,1.0,0);
  col += stars(rd,80,12.0,12,8,1.0,2);
  col += stars(rd,80,12.0,12,8,1.0,3);
  col += stars(rd,80,12.0,12,8,1.0,4);
  col += stars(rd,120,12.0,14,8,1.0,1);
  col += stars(rd,120,12.0,14,8,1.0,2);

  return col;
}



###
###### Ray marching
###

float[10] rayMarch(	vec3 ro, 		// ray origin
                    vec3 rd, 		// ray direction
                    float start,
                    float end) {
  vec3 velocity = rd;
  vec3 p = ro;
  vec4 col = vec4(0);

  for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
    float d = sdfScene(p);
    if (d > end) break;
    float clamp_d = min(d,MAX_DX); //Clamp d to get a continuous looking cloud

    if (sdfDisc(p) < PRECISION) {
      //We hit the disc, we need to go through it and to accumulate color along the way

      vec3 n_vel = normalize(velocity);
      if (p.y > 0) {
        p.y += -3*PRECISION; 
      } else {
        p.y += 3*PRECISION;
      }
      col += disc_color(p);
    } else if (sdfBH(p) < PRECISION) break;

    col += max(clamp_d, 0) * clouds(p); //Adding color from clouds

    p += velocity*d;
  }

  vec3 dir = normalize(velocity);

  float[10] ray_data = float[](
    p.x,
    p.y,
    p.z,
    dir.x,
    dir.y,
    dir.z,
    col.x,
    col.y,
    col.z,
    min(col.w,1)
  );

  return ray_data;
}

void main() {
  vec2 uv = uvCoord;
  
  // We transform the coordinates so they fit the specified aspect ratio and fov
  vec3 uv3d = vec3(uv.x, uv.y, -DIST_TO_CANVAS);
  vec3 col = vec3(0); // Output color vector

  // vec3 rd = normalize(vec3(0.1*uv3d.x, 0.1*uv3d.y/aspect_ratio, uv3d.z)); // ray direction
  vec3 rd = normalize(vec3(fov*uv3d.x, fov*uv3d.y/aspect_ratio, uv3d.z)); // ray direction
  rd = normalize(orientation * rd);

  float[] ray_data = rayMarch(ro, rd, MIN_DIST, MAX_DIST); // distance to sphere
  vec3 p = vec3(ray_data[0],ray_data[1],ray_data[2]);
  vec3 dir = vec3(ray_data[3],ray_data[4],ray_data[5]);
  vec4 added_col = vec4(ray_data[6],ray_data[7],ray_data[8],ray_data[9]);


  if (sdfScene(p) < PRECISION) {
    // The ray hit something
    if (sdfBH(p) < PRECISION) { // The ray hit the black hole
      col = vec3(0);
    }
  } else {
    // Galaxy-like background
    vec3 galaxy_bg = vec3(0.1*pow(cos(dir.y),32));

    // Randomely generated stars through a perlin noise method.
    vec3 rand_stars = create_stars(dir);

    col = rand_stars+galaxy_bg;
  }

  // Output to screen
  float opacity = min(added_col.w,1.0);
  col = opacity*added_col.xyz + (1-opacity)*col;
  fragColor = vec4(col, 1.0);
}
