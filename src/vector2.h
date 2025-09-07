/**
 * @file vector2.h
 * @brief 2D Vector mathematics for ray tracing and geometric calculations
 * 
 * This file provides a simple Vector2 struct with common vector operations
 * needed for ray-surface intersections, reflections, and refractions.
 */

#ifndef VECTOR2_H
#define VECTOR2_H

#include <cmath>

/**
 * @brief 2D Vector structure for geometric calculations
 * 
 * Provides basic vector operations like normalization, length calculation,
 * and supports operator overloading for mathematical operations.
 */
struct Vector2 {
    float x, y;  ///< X and Y components of the vector

    /**
     * @brief Constructor with optional initial values
     * @param x X component (default: 0.0f)
     * @param y Y component (default: 0.0f)
     */
    Vector2(float x = 0.0f, float y = 0.0f) : x(x), y(y) {}

    /**
     * @brief Calculate the magnitude/length of the vector
     * @return The Euclidean length of the vector
     */
    float length() const {
        return std::sqrt(x * x + y * y);
    }

    /**
     * @brief Get a normalized copy of this vector
     * @return A new Vector2 with the same direction but length 1.0
     */
    Vector2 normalized() const {
        float l = length();
        if (l > 0) {
            return Vector2(x / l, y / l);
        }
        return Vector2(0.0f, 0.0f);
    }

    /**
     * @brief Normalize this vector in-place
     * Modifies this vector to have length 1.0 while preserving direction
     */
    void normalize() {
        float l = length();
        if (l > 0) {
            x /= l;
            y /= l;
        }
    }
};

// ============================================================================
// Vector Operations (Operator Overloads and Utility Functions)
// ============================================================================

/**
 * @brief Vector addition
 * @param a First vector
 * @param b Second vector
 * @return Sum of the two vectors
 */
inline Vector2 operator+(const Vector2& a, const Vector2& b) { return Vector2(a.x + b.x, a.y + b.y); }

/**
 * @brief Vector subtraction
 * @param a First vector
 * @param b Second vector
 * @return Difference of the two vectors (a - b)
 */
inline Vector2 operator-(const Vector2& a, const Vector2& b) { return Vector2(a.x - b.x, a.y - b.y); }

/**
 * @brief Scalar multiplication (vector * scalar)
 * @param v Vector to scale
 * @param s Scalar multiplier
 * @return Scaled vector
 */
inline Vector2 operator*(const Vector2& v, float s) { return Vector2(v.x * s, v.y * s); }

/**
 * @brief Scalar multiplication (scalar * vector)
 * @param s Scalar multiplier
 * @param v Vector to scale
 * @return Scaled vector
 */
inline Vector2 operator*(float s, const Vector2& v) { return v * s; }

/**
 * @brief Dot product of two vectors
 * @param a First vector
 * @param b Second vector
 * @return Scalar dot product (a·b)
 */
inline float dot(const Vector2& a, const Vector2& b) { return a.x * b.x + a.y * b.y; }

/**
 * @brief Cross product of two 2D vectors (returns scalar)
 * @param a First vector
 * @param b Second vector
 * @return Scalar cross product (a×b), represents the Z component of 3D cross product
 */
inline float cross(const Vector2& a, const Vector2& b) { return a.x * b.y - a.y * b.x; }

#endif // VECTOR2_H
