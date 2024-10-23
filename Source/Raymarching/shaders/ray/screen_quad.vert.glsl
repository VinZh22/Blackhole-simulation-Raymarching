#version 330 core

// Vertex shader - this code is executed for every vertex of the shape

// Inputs coming from VBOs
layout (location = 0) in vec3 vertex_position; // vertex position in local space (x,y,z)

// Output variables sent to the fragment shader
out vec2 uvCoord;

void main()
{
    uvCoord = vertex_position.xy;

	// gl_Position is a built-in variable which is the expected output of the vertex shader
	gl_Position = vec4(uvCoord, 0.0, 1.0); // gl_Position is the projected vertex position (in normalized device coordinates)
}
