/**
 * @file ray.h
 * @brief Ray class for representing light rays in the simulation
 * 
 * This file defines the Ray class which encapsulates the properties
 * and behavior of light rays, including their origin, direction,
 * color, and rendering functionality.
 */

#ifndef RAY_H
#define RAY_H

#include "vector2.h"

/**
 * @brief Represents a light ray with position, direction and visual properties
 * 
 * The Ray class is used throughout the simulation to represent incident,
 * reflected, refracted, and transmitted light rays. Each ray has an origin
 * point, a direction vector, a color, and a length for visualization.
 */
class Ray {
public:
    Vector2 origin;      ///< Starting point of the ray
    Vector2 direction;   ///< Direction vector (normalized)
    float color[3];      ///< RGB color components (0.0-1.0)
    float length;        ///< Length of the ray for visualization

    /**
     * @brief Constructor for creating a new ray
     * @param origin Starting point of the ray
     * @param direction Direction vector (will be normalized)
     * @param color Array of 3 RGB color components (0.0-1.0)
     * @param length Length of the ray for drawing (default: 2.0f)
     */
    Ray(const Vector2& origin, const Vector2& direction, const float* color, float length = 2.0f);

    /**
     * @brief Render the ray using OpenGL
     * Draws the ray as a colored line from origin to origin + direction * length
     */
    void draw() const;
};

#endif // RAY_H
