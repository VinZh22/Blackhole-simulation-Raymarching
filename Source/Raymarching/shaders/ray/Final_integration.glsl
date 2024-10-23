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

// Ray marching parameters
const int MAX_MARCHING_STEPS = 400;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float PRECISION = 0.001;
const float MAX_DS = 0.1;

//phong reflection const
const float Ka = 0.2;
const float Kd = 0.7;
const float Ks = 0.;
const float specular_exponent = 32.0;

const float BH_RADIUS = .5; // = to 2GM/c² the blackhole's Schwarzschild radius

// Acretion disc parameters
const float DISC_W_IN = 1.85;
const float DISC_W_OUT = 5.0;
const float RING_SPEED = 3.0;

// Black hole parameters

const float PI = 3.14159265359;
const mat3 ID3 = mat3(1,0,0,
                            0,1,0,
                            0,0,1);

//Glass sphere parameters
const float GLASS_RADIUS = 0.5;

//
// Utilities
//

float rand(vec2 co) { //random number generation
    return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

float n (vec3 p) { // Refraction index of space near the blackhole
  return 1/(1-(BH_RADIUS-PRECISION)/length(p));
}

vec3 grad_n (vec3 p) { // Exact gradient of the refraction index
  float n = n(p);
  float r = length(p);
  return -n*n*BH_RADIUS*p/r/r/r;
}

//
//// Scene definition
//

// Primitives

float sdfSphere(vec3 p, vec3 center, float r ) {
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

// Objects definition
float sdfBH(vec3 p) { //The blackhole
  return sdfSphere(p, vec3(0,0,0), 1.0);
}

float sdfDisc (vec3 p) { //The acretion disc
  return max(
    - sdfCappedCylinder(p, 1., DISC_W_IN),
    sdfCappedCylinder(p, 0, DISC_W_OUT)
  );
}

vec3 calcCircleTraj(float r, float vitesse){
  float t = iTime;
  float x = r*cos(vitesse*t);
  float z = r*sin(vitesse*t);
  float y = 0.;
  return vec3(x, y, z);
}

float sdfSceneReflect(vec3 p){
  //give sdfScene of part of the scene that have toal reflection, including the glass sphere (planet)

  float ans = sdfSphere(p, calcCircleTraj(15., 0.5), GLASS_RADIUS); //glass sphere
  ans = min(ans, sdfSphere(p, calcCircleTraj(10., -1.), GLASS_RADIUS-0.2)); //glass sphere
  return ans;
}

float sdfScene(vec3 p) {
  float ans = sdfBH(p);
  ans = min(ans, sdfDisc(p));
  ans = min(ans, sdfSceneReflect(p));
  return ans;
}

int typeScene(vec3 p){
  //return the type of the object that the ray hit
  int ans = 0;
  if(sdfBH(p) < PRECISION){
    ans = 0;
  } else if(sdfDisc(p) < PRECISION){
    ans = 2;
  } else if(sdfSceneReflect(p) < PRECISION){
    ans = 1;
  }
  return ans;
}


//
//// Star background generation
//

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

vec3 stars(
  vec3 rd, //ray direction
  int grid_num, //number of grid points for the stars
  float inverse_width, //The higher this float the smaller the stars
  float scarcity_factor, //The higher this float the more likely a star's intesity is to be reduced
  int inverse_band_width, //The higher this float narrower the generated star band will be
  float global_scaling //A global scaling factor for all the stars intensity
){
  return stars(rd, grid_num, inverse_width, scarcity_factor, inverse_band_width, global_scaling, 1);
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



//
//// Transparent cloud generation
//

vec4 cloud_torus(vec3 p) { // a torus shaped cloud using an sdf
  p = p*vec3(1,2,1);
  float d = sdfTorus(
    p,
    vec2((2*DISC_W_IN+DISC_W_OUT)/3, 0.5),
    vec3(0,0,0),
    ID3
    );
  
  return vec4(0.8,0.8,1,0.03)*max(-d,0);
}

vec4 cloud_sphere(vec3 p) {
  float d = sdfBH(p);
  if (d > 3.0) {return vec4(0);}
  else return vec4(0.8,0.8,1,0.1)*(0.5/(1+10*d*d));
}

vec4 clouds(vec3 p) {
  return cloud_torus(p) + cloud_sphere(p);
}



//
//// Transparent acretion disc generation
//

//We trace a very simple disc and substract lines afterwards
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

  color.w -= lines_neg(r, theta, 10, 5, 1, 5, 5, 2);
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



//
//// Ray marching
//

float[10] rayMarch(	vec3 ro, 		// ray origin
                    vec3 rd, 		// ray direction
                    float end) {
  vec3 velocity = rd;
  vec3 p = ro;
  vec4 col = vec4(0);

  for (int i = 0; i < MAX_MARCHING_STEPS; i++) { //dx = v ds, dv = grad_n(p) * ds, we want dx to be bounded by the sdf distance
    if (length(p-ro) > end) break; //The ray is at infinity, we stop the marching sequence

    float d = sdfScene(p);

    if (sdfDisc(p) < PRECISION) {//The ray hit the disc
      // We make the ray go through the disc
      if (p.y > 0) {
        p.y += -2*PRECISION; 
      } else {
        p.y += 2*PRECISION;
      }
      // We then accumulate color
      col += disc_color(p);
    }

    else if (d < PRECISION) break; //The ray hit something other than the disc, we stop the marching sequence
    
    // Ray propagation according to diffraction laws
    float refractive_idx = n(p);
    float ds = d/length(velocity/refractive_idx);
    ds = min(MAX_DS, ds); //We ensure ds is not too large
    vec3 dx = velocity*ds/refractive_idx;
    p += dx;
    velocity += grad_n(p) * ds;

    col += max(length(dx), 0) * clouds(p); //Adding color from volumetric clouds
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

vec3[2] rayMarch(	vec3 ro, 		// ray origin
                vec3 rd, 		// ray direction
                float start,
                float end) {
  vec3 velocity = rd;
  vec3 p = ro;
  float dt = MAX_DS;

  for (int i = 0; i < MAX_MARCHING_STEPS; i++) { // dx = v dt, dv =  F dt, we want dx to be bounded by the sdf distance, we use variable dt
    float d = sdfScene(p);
    if (d < PRECISION || length(p) > end){
      break;
    }
    
    // Ray propagation according to diffraction laws
    float refractive_idx = n(p);
    float ds = d/length(velocity/refractive_idx);
    ds = min(MAX_DS, ds); //We ensure ds is not too large
    vec3 dx = velocity*ds/refractive_idx;
    p += dx;
    velocity += grad_n(p) * ds;
  }

  vec3[2] pos_dir;
  pos_dir[0] = p;
  pos_dir[1] = normalize(velocity);

  return pos_dir;
}

//
//// Reflections
//

vec3 phong(vec3 normal, vec3 lightDirection, vec3 cameraDirection, vec3 color_object, vec3 color_light) {
  // Phong coefficient (diffuse, specular)
  // Copied and adapted the code from what we did on TD2 //

  // Unit direction toward the light
  vec3 L = normalize(lightDirection);
  vec3 N = normalize(normal);

  // Diffuse coefficient
  float diffuse_component = max(dot(N,L),0.0);

  // Specular coefficient
  float specular_component = 0.0;
  if(diffuse_component>0.0){
    vec3 R = reflect(-L,N); // reflection of light vector relative to the normal.
    vec3 V = normalize(cameraDirection);
    specular_component = pow( max(dot(R,V),0.0), specular_exponent );
  }
  vec3 color_shading = (Ka + Kd * diffuse_component) * color_object + Ks * specular_component * color_light;

  return color_shading;
}

vec3 estimateNormal(vec3 p) {
    return normalize(vec3(
        sdfScene(vec3(p.x + PRECISION, p.y, p.z)) - sdfScene(vec3(p.x - PRECISION, p.y, p.z)),
        sdfScene(vec3(p.x, p.y + PRECISION, p.z)) - sdfScene(vec3(p.x, p.y - PRECISION, p.z)),
        sdfScene(vec3(p.x, p.y, p.z  + PRECISION)) - sdfScene(vec3(p.x, p.y, p.z - PRECISION))
    ));
}


//
//// Main function
//

void main() {
  vec2 uv = uvCoord;
  
  // We transform the coordinates so they fit the specified aspect ratio and fov
  vec3 uv3d = vec3(uv.x, uv.y, -DIST_TO_CANVAS);
  vec3 col = vec3(1,0,0); // Output color vector
  vec3 rd = normalize(vec3(fov*uv3d.x, fov*uv3d.y/aspect_ratio, uv3d.z)); // ray direction
  rd = normalize(orientation * rd);

  float[10] ray_data = rayMarch(ro, rd, MAX_DIST); // results of ray marching operation
  vec3 p = vec3(ray_data[0],ray_data[1],ray_data[2]);
  vec3 dir = vec3(ray_data[3],ray_data[4],ray_data[5]);
  vec4 added_col = vec4(ray_data[6],ray_data[7],ray_data[8],ray_data[9]);


  if (sdfScene(p) < PRECISION) {
    // The ray hit something
    int type = typeScene(p);
    if (type == 0) { // The ray hit the black hole
      col = vec3(0);
    }
    else if (type == 1) { // The ray hit the glass sphere
      vec3 normal = estimateNormal(p);
      vec3 lightDirection = normalize(dir); // direction to the light source
      vec3 cameraDirection = normalize(ro - p); // direction to the camera

      vec3 reflVec = reflect(lightDirection, normal);
      // we do another raymarch from the reflection point

      vec3[2] pos_dir_refl = rayMarch(p+reflVec*MAX_DS, reflVec, MIN_DIST, MAX_DIST); // position and direction of ray at end of raymarch
      float d_refl = sdfScene(pos_dir_refl[0]); // distance to sdf
      vec3 objectColor = vec3(0.); vec3 lightColor = vec3(1.);
      if (d_refl > PRECISION) {
        // The reflected ray didn't hit anything: We display stars at infinity.

        // Galaxy-like background
        vec3 galaxy_bg = vec3(0.1*pow(cos(pos_dir_refl[1].y),32));

        // Randomely generated stars through a perlin noise method.
        vec3 rand_stars = create_stars(pos_dir_refl[1]);
        objectColor = rand_stars+galaxy_bg;

      } else {
        //it touched something
          int typeRef = typeScene(pos_dir_refl[0]);
          if (typeRef == 0) { // The ray hit the black hole
            objectColor = vec3(0);
          }
          else if (typeRef == 2){ // The ray hit the acretion disc
            objectColor = disc_color(pos_dir_refl[0]).xyz;
          }
          else if (typeRef==1){
            objectColor = vec3(1);
          }
      }
      lightDirection = normalize(ro-p); // triche : direction to the light source
      col = phong(normal, lightDirection, cameraDirection, objectColor, lightColor);
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