/**
 * @file ray.cpp
 * @brief Implementation of the Ray class for light ray visualization
 */

#include "ray.h"
#include "simulation.h"

/**
 * @brief Constructor implementation
 * Creates a new ray with specified properties and normalizes the direction vector
 */
Ray::Ray(const Vector2& origin, const Vector2& direction, const float* color, float length)
    : origin(origin), direction(direction.normalized()), length(length) {
    // Copy RGB color components
    this->color[0] = color[0];
    this->color[1] = color[1];
    this->color[2] = color[2];
}

/**
 * @brief Draw the ray as a colored line
 * Renders the ray from its origin to the endpoint (origin + direction * length)
 */
void Ray::draw() const {
    Vector2 end = origin + direction * length;
    draw_line(origin.x, origin.y, end.x, end.y, color[0], color[1], color[2], 5.0f);
}
