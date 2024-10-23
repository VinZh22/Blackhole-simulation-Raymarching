#include "ray_source.hpp"
#include "environment.hpp"

using namespace cgp;

mat3 ray_source::orientation_mat() {
    return orientation().matrix();
}

void ray_source::send_uniform (environment_structure& environment) {
    environment.uniform_generic.uniform_vec3["ro"] = position();
    environment.uniform_generic.uniform_mat3["orientation"] = orientation_mat();
    environment.uniform_generic.uniform_float["aspect_ratio"] = aspect_ratio;
    environment.uniform_generic.uniform_float["fov"] = fov;
    environment.uniform_generic.uniform_float["DIST_TO_CANVAS"] = DIST_TO_CANVAS;
}

void ray_source_controller::action_mouse_move()
	{
		// Preconditions
		assert_cgp_no_msg(inputs != nullptr);
		assert_cgp_no_msg(window != nullptr);
		if (!is_active) return;

		vec2 const& p1 = inputs->mouse.position.current;
		vec2 const& p0 = inputs->mouse.position.previous;
		vec2 const dp = p1 - p0;

		bool const event_valid = !inputs->mouse.on_gui;
		bool const click_left = inputs->mouse.click.left;
		bool const click_right = inputs->mouse.click.right;
		bool const ctrl = inputs->keyboard.ctrl;

		if (event_valid) { // If the mouse cursor is not on the ImGui area

			if (click_left && !ctrl)     // Rotation of the camera around its center
				camera_model.manipulator_rotate_roll_pitch_yaw(0, dp.y, -dp.x);
			else if (click_left && ctrl) // Translate/Pan the camera in the viewspace plane
				camera_model.manipulator_translate_in_plane(p1 - p0);
			else if (click_right && !ctrl) // Move the camera closer/further with respect to its center
				camera_model.manipulator_scale_distance_to_center((p1 - p0).y);
			else if (click_right && ctrl) // Translate the camera center in front/back
				camera_model.manipulator_translate_front((p1 - p0).y);
		}
	}