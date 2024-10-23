#version 330 core

//
////// Ray marching shader
// Implements star generation with perlin noise and bended light rays (not relativistic)

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
const float MIN_DT = 0.1;
const float PRECISION = 0.001;

const float PI = 3.14159265359;

//phong reflection const
const float Ka = 0.2;
const float Kd = 0.7;
const float Ks = 0.;
const float specular_exponent = 32.0;
//blackhole
const float BH_RADIUS = 0.1; // = to c²/2GM the blackhole radius
const float BH_FORCE = 0.04;
const vec3 BH_CENTER = vec3(0.);
//glass sphere
const float GLASS_RADIUS = 1.;
const vec3 GLASS_CENTER = vec3(10., 0., -2.);


float rand(vec2 co) {
    return fract(sin(dot(co.xy ,vec2(12.9898,78.233))) * 43758.5453);
}

float sdSphere(vec3 p, vec3 center, float r ) {
  return length(p-center) - r;
}

float sdSphere(vec3 p, float r ) {
  //overload considering the center is 0
  return length(p) - r;
}

float sdSceneBH(vec3 p){
  //give sdScene of the blackhole
  return sdSphere(p, BH_CENTER, BH_RADIUS); //blackhole;
}

vec3 calcCircleTraj(float r, float vitesse){
  float t = iTime;
  float x = r*cos(vitesse*t);
  float z = r*sin(vitesse*t);
  float y = 0.;
  return vec3(x, y, z);
}

float sdSceneReflect(vec3 p){
  //give sdScene of part of the scene that have toal reflection, including the glass sphere (planet)

  float ans = sdSphere(p, calcCircleTraj(5., 1.), GLASS_RADIUS); //glass sphere
  ans = min(ans, sdSphere(p, calcCircleTraj(10., 4.), GLASS_RADIUS-0.5)); //glass sphere
  return ans;
}

float sdSceneReste(vec3 p){
  float ans = MAX_DIST;//ship(p);
  return ans;
}

float sdScene(vec3 p) {
  float ans = sdSceneBH(p);
  ans = min(ans, sdSceneReste(p));
  ans = min(ans, sdSceneReflect(p));
  return ans;
}


int typeScene(vec3 p){
  // 0 = blackhole, 1 = glass sphere, 2 = reste
  int ans = 0;
  float val = sdSceneBH(p);
  if (val > sdSceneReflect(p)){
    ans = 1;
    val = sdSceneReflect(p);
  }
  if (val > sdSceneReste(p)){
    ans = 2;
    val = sdSceneReste(p);
  }
  return ans;
}

vec3 force (vec3 p) {
  float r = length(p-BH_CENTER);
  return -normalize(p) * BH_FORCE / (r*r);
}

vec3[2] rayMarch(	vec3 ro, 		// ray origin
                vec3 rd, 		// ray direction
                float start,
                float end) {
  vec3 velocity = rd;
  vec3 p = ro;
  float dt = MIN_DT;

  for (int i = 0; i < MAX_MARCHING_STEPS; i++) { // dx = v dt, dv =  F dt, we want dx to be bounded by the sdf distance, we use variable dt
    float d = sdScene(p);
    if (d < PRECISION || length(p) > end){
      break;
    }
    
    dt = d/length(velocity);
    dt = max(MIN_DT, dt);

    p += velocity*dt;
    velocity += force(p) * dt;
  }

  vec3[2] pos_dir;
  pos_dir[0] = p;
  pos_dir[1] = normalize(velocity);

  return pos_dir;
}

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
        sdScene(vec3(p.x + PRECISION, p.y, p.z)) - sdScene(vec3(p.x - PRECISION, p.y, p.z)),
        sdScene(vec3(p.x, p.y + PRECISION, p.z)) - sdScene(vec3(p.x, p.y - PRECISION, p.z)),
        sdScene(vec3(p.x, p.y, p.z  + PRECISION)) - sdScene(vec3(p.x, p.y, p.z - PRECISION))
    ));
}


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
  vec3 col = stars(rd,10,9.0,8.0,4,2.0, 1);
  col += stars(rd,20,9.0,10,6,1.0, 1);
  col += stars(rd,20,9.0,10,6,1.0,2);
  col += stars(rd,40,12.0,12,6,1.0, 1);
  col += stars(rd,40,12.0,12,6,1.0,2);
  col += stars(rd,40,12.0,12,6,1.0,3);
  col += stars(rd,40,12.0,12,6,1.0,4);
  col += stars(rd,80,12.0,12,8,1.0, 1);
  col += stars(rd,80,12.0,12,8,1.0,2);
  col += stars(rd,80,12.0,12,8,1.0,3);
  col += stars(rd,80,12.0,12,8,1.0,4);
  col += stars(rd,120,12.0,14,8,1.0, 1);
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
    int i = typeScene(pos_dir[0]);
    if (i==0){
      //black hole
      col = vec3(0.);
    }
    else if (i==1){
      //glass sphere
      vec3 normal = estimateNormal(pos_dir[0]); // normal at the hit
      vec3 lightDirection = normalize(pos_dir[1]); // direction to the light source
      vec3 cameraDirection = normalize(ro - pos_dir[0]); // direction to the camera

      vec3 reflVec = reflect(pos_dir[1], normal);
      // we do another raymarch from the reflection point

      vec3[2] pos_dir_refl = rayMarch(pos_dir[0]+reflVec*MIN_DT, reflVec, MIN_DIST, MAX_DIST); // position and direction of ray at end of raymarch
      float d_refl = sdScene(pos_dir_refl[0]); // distance to sdf
      vec3 objectColor = vec3(0.); vec3 lightColor = vec3(1.);
      if (d_refl > PRECISION) {
        // The reflected ray didn't hit anything: We display stars at infinity.

        // Galaxy-like background
        vec3 galaxy_bg = vec3(0.1*pow(cos(pos_dir_refl[1].y),32));

        // Randomely generated stars through a perlin noise method.
        vec3 rand_stars = create_stars(pos_dir_refl[1]);
        objectColor = rand_stars+galaxy_bg;

      } else {
        //it touched the black hole
        objectColor = vec3(0.);
      }
      lightDirection = normalize(ro-pos_dir[0]); // direction to the light source
      col = phong(normal, lightDirection, cameraDirection, objectColor, lightColor);
    }
    else{
      //reste
      vec3 normal = estimateNormal(pos_dir[0]); // normal at the hit
      vec3 lightDirection = normalize(ro - pos_dir[0]); // direction to the light source
      vec3 cameraDirection = normalize(ro - pos_dir[0]); // direction to the camera
      //tout est en rouge pour l'instant (ya rien de tt facon)
      vec3 objectColor = vec3(1., 0., 0.); vec3 lightColor = vec3(1.);
      col = phong(normal, lightDirection, cameraDirection, objectColor, lightColor);
    }
  }

  // Output to screen
  fragColor = vec4(col, 1.0);
}
