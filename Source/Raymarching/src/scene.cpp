#include "scene.hpp"
#include "ray_source.hpp"

using namespace cgp;



// This function is called only once at the beginning of the program
// This function can contain any complex operation that can be pre-computed once
void scene_structure::initialize()
{
	std::cout << "Start function scene_structure::initialize()" << std::endl;

	// Set the behavior of the camera and its initial position
	// ********************************************** //
	ray_source.initialize(inputs, window); 
	ray_source.set_rotation_axis_z(); // camera rotates around z-axis
	//   look_at(camera_position, targeted_point, up_direction)
	ray_source.look_at(
		{ 0.0f, 1.0f, 3.001f } /* position of the camera in the 3D scene */,
		{0,0,3.0f} /* targeted point in 3D scene */,
		{0,1,0} /* direction of the "up" vector */);

	ray_source.camera_model.center_of_rotation = {0.0f, 1.0f, 3.0f};
	ray_source.camera_model.distance_to_center = 0.001f;

	// General information
	display_info();

	// Create the shapes seen in the 3D scene and initialize the on the gpu
	// ********************************************** //
	display_canvas.initialize_data_on_gpu(mesh_primitive_quadrangle({-1.0f,-1.0f,0},{1.0f,-1.0f,0},{1.0f,1.0f,0},{-1.0f,1.0f,0}));


	// Load a shader from a file
	opengl_shader_structure custom_shader;
	custom_shader.load(
		project::path + "shaders/ray/screen_quad.vert.glsl", 
		project::path + "shaders/ray/Final_integration.glsl");

	// Affect the loaded shader to the mesh_drawable
	display_canvas.shader = custom_shader;
	
	// TODO code a way to create a cubemap texture and send it to the shader as a uniform array;

	std::cout << "End function scene_structure::initialize()" << std::endl;

}


// This function is called permanently at every new frame
// Note that you should avoid having costly computation and large allocation defined there. This function is mostly used to call the draw() functions on pre-existing data.
void scene_structure::display_frame()
{
	// Get the current aspect ratio of the window
	ray_source.camera_model.aspect_ratio = window.aspect_ratio();
	ray_source.camera_model.send_uniform(environment);

	// Set the light to the current position of the camera
	environment.light = ray_source.camera_model.position();

	// Update time
	timer.update();
	environment.uniform_generic.uniform_float["iTime"] = timer.t;
	

	// the general syntax to display a mesh is:
	//   draw(mesh_drawableName, environment);
	// We display the canvas used to perform ray tracing
	draw(display_canvas, environment);
}

void scene_structure::display_gui()
{

}

void scene_structure::display_info()
{
	std::cout << "\nCAMERA CONTROL:" << std::endl;
	std::cout << "-----------------------------------------------" << std::endl;
	std::cout << ray_source.doc_usage() << std::endl;
	std::cout << "-----------------------------------------------\n" << std::endl;


	std::cout << "\nSCENE INFO:" << std::endl;
	std::cout << "-----------------------------------------------" << std::endl;
	std::cout << "Introductory scene for INF443." << std::endl;
	std::cout << "-----------------------------------------------\n" << std::endl;
}

void scene_structure::mouse_move_event()
{
	if (!inputs.keyboard.shift)
		ray_source.action_mouse_move();
}
void scene_structure::mouse_click_event()
{
	ray_source.action_mouse_click(environment.dummy_camera_view);
}
void scene_structure::keyboard_event()
{
	ray_source.action_keyboard(environment.dummy_camera_view);
}
void scene_structure::idle_frame()
{
	ray_source.idle_frame(environment.dummy_camera_view);
}