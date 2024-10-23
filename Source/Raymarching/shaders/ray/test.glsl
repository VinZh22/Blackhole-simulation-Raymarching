#version 330 core

//
////// Ray-Tracing shader to be implemented
// 

// Inputs of th ray shader - fragment position
in vec2 uvCoord;

// Output of the ray shader - fragment color
out vec4 fragColor;

// Uniform values sent by the CGP based c++ code
uniform float time;     // time since simulation start

uniform float aspect_ratio;

vec3 sdfSphere(vec3 p, float r) {
  float d = length(p) - r;

  return d > 0. ? vec3(1.0) : vec3(0., 0., 1.);
}

//Todo : implement ways to account for aspect ratio.

void main()
{
  vec2 uv = uvCoord;
  // uv = uv / vec2(1.0, aspect_ratio); // Fix aspect ratio
  vec3 uv3d = vec3(uv.x, uv.y/aspect_ratio, 0.0);

  vec3 col = sdfSphere(uv3d, 0.2); // Call this function on each pixel to check if the coordinate lies inside or outside of the circle

  // Output to screen
  fragColor = vec4(col,1.0);
}