/**
 * @file refraction.cpp
 * @brief Implementation of light refraction simulation
 * 
 * This file handles the visualization of light bending when passing between
 * different media (e.g., air to water), demonstrating Snell's law of refraction
 * and total internal reflection.
 */

#include "simulation.h"
#include "optics.h"
#include "ray.h"
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
    const OpticalSurface interfaceSurface = makeFlatSurface(
        SurfaceType::FlatInterface,
        Vector2(0.0f, interfaceY),
        Vector2(1.0f, 0.0f),
        Vector2(0.0f, 1.0f),
        -1.0f,
        1.0f,
        makeDielectricMaterial(refractiveIndex1),
        makeDielectricMaterial(refractiveIndex2)
    );
    
    // Draw interface line between two media (horizontal gray line)
    draw_line(-1.0f, interfaceY, 1.0f, interfaceY, 0.7f, 0.7f, 0.7f, 3.0f);

    // Convert current rotation angle to radians
    float rad = rotationAngle * M_PI / 180.0f;
    Vector2 inc_dir(-std::cos(rad), -std::sin(rad));

    // Calculate starting point of incident ray (offset from interface)
    Vector2 target_point(rayPosition, interfaceY);
    Vector2 start_point = target_point - inc_dir * 1.5f;
    SurfaceHit hit = intersectRayWithSurface(start_point, inc_dir, interfaceSurface);
    if (!hit.hit) {
        return;
    }

    // Draw incident ray approaching the interface
    Ray incident_ray(start_point, inc_dir, rayAColor, hit.distance);
    incident_ray.draw();

    // Draw normal line (dashed green vertical line) at intersection point
    draw_dashed_line(hit.point.x, interfaceY - 0.5f, hit.point.x, interfaceY + 0.5f,
                     0.0f, 1.0f, 0.0f, 1.0f);

    Vector2 reflected_dir = reflectDirection(inc_dir, hit.normal);
    RefractionResult refraction = refractDirection(
        inc_dir,
        hit.normal,
        hit.incomingMaterial.refractiveIndex,
        hit.outgoingMaterial.refractiveIndex
    );
    totalInternalReflection = refraction.totalInternalReflection;
    if (totalInternalReflection) {
        float reflectedColor[] = {
            rayCColor[0],
            rayCColor[1],
            rayCColor[2]
        };
        Ray reflected_ray(hit.point, reflected_dir, reflectedColor, 2.0f);
        reflected_ray.draw();
    } else {
        float reflectionBrightness = 0.4f;
        float reflectedColor[] = {
            rayCColor[0] * reflectionBrightness,
            rayCColor[1] * reflectionBrightness,
            rayCColor[2] * reflectionBrightness
        };
        Ray reflected_ray(hit.point, reflected_dir, reflectedColor, 2.0f);
        reflected_ray.draw();
        Ray refracted_ray(hit.point, refraction.direction, rayBColor, 2.0f);
        refracted_ray.draw();
    }

    // Draw intersection point as a white dot
    draw_point(hit.point.x, hit.point.y, 1.0f, 1.0f, 1.0f, 6.0f);
}
