/**
 * @file refraction.cpp
 * @brief Implementation of light refraction simulation
 * 
 * This file handles the visualization of light bending when passing between
 * different media (e.g., air to water), demonstrating Snell's law of refraction
 * and total internal reflection.
 */

#include "simulation.h"
#include "ray.h"
#include "vector2.h"
#include <cmath>

/**
 * @brief Draw the refraction simulation scene
 * 
 * This function renders:
 * - A horizontal interface between two media
 * - An incident light ray approaching the interface
 * - A refracted ray bending according to Snell's law
 * - A reflected ray (always present, varies in intensity)
 * - Normal line perpendicular to the interface
 * - Handles total internal reflection when applicable
 * 
 * Uses Snell's law: n1 * sin(θ1) = n2 * sin(θ2)
 */
void draw_refraction() {
    const float interfaceY = 0.0f;  // Y coordinate of the media interface
    
    // Draw interface line between two media (horizontal gray line)
    draw_line(-1.0f, interfaceY, 1.0f, interfaceY, 0.7f, 0.7f, 0.7f, 3.0f);

    // Convert current rotation angle to radians
    float rad = rotationAngle * M_PI / 180.0f;
    Vector2 inc_dir(-std::cos(rad), -std::sin(rad));

    // Calculate starting point of incident ray (offset from interface)
    Vector2 start_point(rayPosition, interfaceY);
    start_point = start_point - inc_dir * 1.5f;

    // Draw incident ray approaching the interface
    Ray incident_ray(start_point, inc_dir, rayAColor, 1.5f);
    incident_ray.draw();

    // Draw normal line (dashed green vertical line) at intersection point
    draw_dashed_line(rayPosition, interfaceY - 0.5f, rayPosition, interfaceY + 0.5f,
                     0.0f, 1.0f, 0.0f, 1.0f);

    // Calculate incident angle from the normal
    float incident = std::fabs(90.0f - std::fmod(rotationAngle, 180.0f));
    float incidentRad = incident * M_PI / 180.0f;
    
    // Apply Snell's law: n1 * sin(θ1) = n2 * sin(θ2)
    // Solve for sin(θ2) = (n1/n2) * sin(θ1)
    float sinRefr = (n1 / n2) * sinf(incidentRad);

    // Check for total internal reflection
    if (sinRefr > 1.0f || sinRefr < -1.0f) {
        // Total internal reflection occurs - no transmitted ray
        totalInternalReflection = true;
    } else {
        // Normal refraction occurs
        totalInternalReflection = false;

        // Calculate and draw reflected ray (always present)
        Vector2 refl_dir(inc_dir.x, -inc_dir.y);  // Mirror Y component for reflection
        float reflectionBrightness = totalInternalReflection ? 1.0f : 0.4f;
        float reflectedColor[] = { 
            rayCColor[0] * reflectionBrightness, 
            rayCColor[1] * reflectionBrightness, 
            rayCColor[2] * reflectionBrightness 
        };
        Ray reflected_ray(Vector2(rayPosition, interfaceY), refl_dir, reflectedColor, 2.0f);
        reflected_ray.draw();

        // Calculate refracted ray direction using Snell's law
        float refractionAngle = asinf(sinRefr);
        float refX = copysignf(sinf(refractionAngle), inc_dir.x);  // Preserve horizontal direction
        float refY = copysignf(cosf(refractionAngle), inc_dir.y);  // Bend according to Snell's law
        Vector2 refr_dir = { refX, refY };
        refr_dir.normalize();
        
        // Draw refracted ray bending into the second medium
        Ray refracted_ray(Vector2(rayPosition, interfaceY), refr_dir, rayBColor, 2.0f);
        refracted_ray.draw();
    }

    // Draw intersection point as a white dot
    draw_point(rayPosition, interfaceY, 1.0f, 1.0f, 1.0f, 6.0f);
}