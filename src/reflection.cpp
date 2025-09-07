/**
 * @file reflection.cpp
 * @brief Implementation of mirror reflection simulation
 * 
 * This file handles the visualization of light reflection off a mirror surface,
 * demonstrating the law of reflection: angle of incidence equals angle of reflection.
 */

#include "simulation.h"
#include "ray.h"
#include "vector2.h"
#include <cmath>

/**
 * @brief Draw the reflection simulation scene
 * 
 * This function renders:
 * - A horizontal mirror surface
 * - An incident light ray approaching the mirror
 * - A reflected light ray bouncing off the mirror
 * - Normal line (dashed) perpendicular to the mirror
 * - Intersection point where the ray hits the mirror
 * 
 * The reflection follows the law: angle of incidence = angle of reflection
 */
void draw_reflection() {
    const float mirrorY = -0.8f;  // Y coordinate of the mirror surface

    // Draw the mirror as a horizontal gray line
    draw_line(-1.0f, mirrorY, 1.0f, mirrorY, 0.7f, 0.7f, 0.7f, 3.0f);

    // Constrain rotation angle to 0-360 degrees
    float constrainedAngle = fmod(rotationAngle, 360.0f);
    if (constrainedAngle < 0) constrainedAngle += 360.0f;

    // Convert angle to radians and calculate incident direction vector
    float rad = constrainedAngle * M_PI / 180.0f;
    Vector2 incident_dir(-std::cos(rad), -std::sin(rad));
    
    // Calculate intersection point on the mirror surface
    Vector2 intersection_point(rayPosition, mirrorY);
    
    // Calculate starting point of incident ray (offset backwards from intersection)
    Vector2 start_point = intersection_point - incident_dir * 2.0f;

    // Draw incident ray (Ray A) - approaching the mirror
    Ray incident_ray(start_point, incident_dir, rayAColor, 2.0f);
    incident_ray.draw();

    // Draw normal line (dashed green vertical line) at intersection point
    draw_dashed_line(intersection_point.x, mirrorY, intersection_point.x, 1.0f,
                     0.0f, 1.0f, 0.0f, 1.0f);

    // Calculate reflected direction (mirror the Y component, keep X component)
    Vector2 reflected_dir(incident_dir.x, -incident_dir.y);

    // Draw reflected ray (Ray B) - bouncing off the mirror
    Ray reflected_ray(intersection_point, reflected_dir, rayBColor, 2.0f);
    reflected_ray.draw();

    // Draw intersection point as a white dot
    draw_point(intersection_point.x, intersection_point.y, 1.0f, 1.0f, 1.0f, 6.0f);
}