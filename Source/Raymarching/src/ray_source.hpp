#pragma once

#include "environment.hpp"
#include "cgp/cgp.hpp"

using namespace cgp;

// Structure used to define the source of the rays used in ray tracing based shaders
struct ray_source : public camera_orbit_euler {
    float aspect_ratio = 1.0;
    float fov = 1.0;

    void send_uniform(environment_structure& environment); // Sends uniform parameters with ray origin informations to open gl
    
    private:
        float width;
        float height;
        const float DIST_TO_CANVAS = 1.0f;
    mat3 orientation_mat(); // Returns the ray direction of the camera
};

struct ray_source_controller : public camera_controller_orbit_euler {
    ray_source camera_model;

    //TODO add functions to move camera around using mouse and keys drawing inspiration from camera_controller_orbit ?
    void action_mouse_move();
};