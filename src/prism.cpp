/**
 * @file prism.cpp
 * @brief Implementation of complex prism light simulation
 * 
 * This file handles the most complex simulation mode, showing how light interacts
 * with various prism shapes including rectangles, triangles, circles, and semicircles.
 * It demonstrates multiple phenomena:
 * - Ray-surface intersections with geometric shapes
 * - Multiple internal reflections and refractions
 * - Light dispersion (basic RGB separation)
 * - Fresnel reflection coefficients
 * - Total internal reflection in complex geometries
 */

#include "simulation.h"
#include "vector2.h"
#include "ray.h"
#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <vector>

// ============================================================================
// Global Variables and Constants
// ============================================================================

// Wavelengths for different colors (illustrative, not currently used for dispersion)
static float wavelengths[3] = { 700e-9f, 550e-9f, 450e-9f }; // Red, Green, Blue
extern float refractiveIndex2_RGB[3];  // RGB-specific refractive indices

// Prism type enumeration (local copy for this file)
enum PrismType { RECTANGLE, TRIANGLE, CIRCLE, SEMICIRCLE };

// Prism properties
PrismType currentPrismType = RECTANGLE;    ///< Current prism shape
float prismRotation = 0.0f;                ///< Prism rotation angle in degrees
float prismRectHeight = 0.3f;              ///< Height of rectangular prism

// External references to global simulation variables
extern float rotationAngle;                 ///< Incident ray angle
extern float refractiveIndex1;              ///< Environmental medium refractive index
extern float refractiveIndex2;              ///< Prism material refractive index
extern float rayAColor[3];                  ///< Incident ray color (red channel)
extern float rayBColor[3];                  ///< Internal/transmitted ray color (green channel)
extern float rayCColor[3];                  ///< Reflected ray color (blue channel)

// Ray tracking and display variables
float rayIncidentAngle = 0.0f;             ///< Incident angle for UI display
float rayInternalAngle = 0.0f;             ///< Internal ray angle for UI display
float rayExitAngle = 0.0f;                 ///< Exit ray angle for UI display
float rayDistance = 0.0f;                  ///< Total ray travel distance
std::vector<std::string> allRayAngles;     ///< Log of all ray interaction events
int internalRayCount = 0;                  ///< Count of internal ray bounces

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Clamp a value between minimum and maximum bounds
 * @param value Input value to clamp
 * @param min Minimum allowed value
 * @param max Maximum allowed value
 * @return Clamped value
 */
inline float clamp(float value, float min, float max) {
    return (value < min) ? min : ((value > max) ? max : value);
}

/**
 * @brief Calculate Fresnel reflection coefficient
 * Determines how much light is reflected vs transmitted at an interface
 * @param cosI Cosine of incident angle
 * @param n1 Refractive index of first medium
 * @param n2 Refractive index of second medium
 * @return Reflection coefficient (0.0 = no reflection, 1.0 = total reflection)
 */
float fresnel(float cosI, float n1, float n2) {
    float sinT2 = (n1/n2)*(n1/n2)*(1 - cosI*cosI);
    if (sinT2 > 1.0f) return 1.0f; // Total internal reflection
    
    float cosT = std::sqrt(1 - sinT2);
    // Calculate s-polarized and p-polarized reflection coefficients
    float rs = (n1*cosI - n2*cosT) / (n1*cosI + n2*cosT);
    float rp = (n2*cosI - n1*cosT) / (n2*cosI + n1*cosT);
    
    // Return average of both polarizations
    return 0.5f*(rs*rs + rp*rp);
}

// Forward declarations for drawing functions
void draw_line(float x1, float y1, float x2, float y2, float r, float g, float b, float width);
void draw_point(float x, float y, float r, float g, float b, float size);

/**
 * @brief Check if a ray intersects with a line segment
 * Used for ray-prism surface intersection calculations
 * @param rayOrigin Starting point of the ray
 * @param rayDir Direction vector of the ray (should be normalized)
 * @param segStart Start point of the line segment
 * @param segEnd End point of the line segment
 * @param intersection Output: intersection point if found
 * @param t Output: parameter along ray where intersection occurs
 * @return true if intersection exists, false otherwise
 */
bool raySegmentIntersection(const Vector2& rayOrigin, const Vector2& rayDir, 
                           const Vector2& segStart, const Vector2& segEnd, 
                           Vector2& intersection, float& t) {
    Vector2 segDir = segEnd - segStart;
    float crossProduct = cross(rayDir, segDir);
    
    // Check if ray and segment are parallel (no intersection)
    if (std::abs(crossProduct) < 1e-6f)
        return false;
    
    Vector2 delta = segStart - rayOrigin;
    float t1 = cross(delta, segDir) / crossProduct;
    float t2 = cross(delta, rayDir) / crossProduct;
    
    // Check if intersection occurs within the ray (t1 >= 0) and segment (0 <= t2 <= 1)
    if (t1 >= 0.0f && t2 >= 0.0f && t2 <= 1.0f) {
        t = t1;
        intersection = rayOrigin + rayDir * t1;
        return true;
    }
    
    return false;
}

/**
 * @brief Draw a dashed line for visual guides
 * @param x0 Start X coordinate
 * @param y0 Start Y coordinate  
 * @param dx Direction X component (normalized)
 * @param dy Direction Y component (normalized)
 * @param length Total length of the dashed line
 * @param r Red color component
 * @param g Green color component
 * @param b Blue color component
 * @param width Line width
 */
void draw_dashed_line(float x0, float y0, float dx, float dy, float length, 
                     float r, float g, float b, float width) {
    int dashCount = 16;
    float dashLength = length / dashCount;
    
    // Draw every other segment to create dashed effect
    for (int i = 0; i < dashCount; i += 2) {
        float t0 = i * dashLength;
        float t1 = (i + 1) * dashLength;
        float startX = x0 + dx * t0;
        float startY = y0 + dy * t0;
        float endX = x0 + dx * t1;
        float endY = y0 + dy * t1;
        draw_line(startX, startY, endX, endY, r, g, b, width);
    }
}

/**
 * @brief Draw the outline of the current prism shape
 * Renders the prism boundary based on currentPrismType and applies rotation
 */
void draw_prism_shape() {
    glPushMatrix();
    glTranslatef(0.0f, 0.0f, 0.0f);
    glRotatef(prismRotation, 0.0f, 0.0f, 1.0f);  // Apply prism rotation
    glColor3f(0.8f, 0.8f, 0.8f);                 // Gray outline color
    
    
    switch (currentPrismType) {
        case RECTANGLE: {
            // Draw rectangular prism outline
            glBegin(GL_LINE_LOOP);
            glVertex2f(-0.3f, -prismRectHeight);
            glVertex2f(0.3f, -prismRectHeight);
            glVertex2f(0.3f, prismRectHeight);
            glVertex2f(-0.3f, prismRectHeight);
            glEnd();
            break;
        }
        case TRIANGLE: {
            // Draw triangular prism outline (equilateral triangle)
            glBegin(GL_LINE_LOOP);
            glVertex2f(0.0f, 0.3f);      // Top vertex
            glVertex2f(0.3f, -0.3f);     // Bottom right vertex
            glVertex2f(-0.3f, -0.3f);    // Bottom left vertex
            glEnd();
            break;
        }
        case CIRCLE: {
            // Draw circular prism outline
            glBegin(GL_LINE_LOOP);
            for (int i = 0; i < 100; ++i) {
                float angle = i * 2.0f * M_PI / 100.0f;
                glVertex2f(0.3f * std::cos(angle), 0.3f * std::sin(angle));
            }
            glEnd();
            break;
        }
        case SEMICIRCLE: {
            // Draw semicircular prism outline
            // Draw the curved top part
            glBegin(GL_LINE_STRIP);
            for (int i = 0; i <= 50; ++i) {
                float angle = M_PI * i / 50.0f;  // 0 to π radians
                glVertex2f(0.3f * std::cos(angle), 0.3f * std::sin(angle));
            }
            glEnd();
            
            // Draw the flat bottom border
            glBegin(GL_LINES);
            glVertex2f(-0.3f, 0.0f);
            glVertex2f(0.3f, 0.0f);
            glEnd();
            break;
        }
    }
    glPopMatrix();
}

// ============================================================================
// Physics and Geometry Functions
// ============================================================================

/**
 * @brief Calculate reflection vector from incident ray and surface normal
 * Implements the law of reflection: angle of incidence = angle of reflection
 * @param incident Incident ray direction vector
 * @param normal Surface normal vector (should point towards the incident medium)
 * @return Reflected ray direction vector
 */
Vector2 reflect(const Vector2& incident, const Vector2& normal) {
    float dotProduct = dot(incident, normal);
    return incident - 2.0f * dotProduct * normal;
}

/**
 * @brief Calculate refracted ray direction using Snell's law
 * Handles both normal refraction and total internal reflection
 * @param incident Incident ray direction vector
 * @param normal Surface normal vector (pointing towards incident medium)
 * @param iorRatio Ratio of refractive indices (n1/n2)
 * @param totalInternalReflection Output: set to true if TIR occurs
 * @return Refracted ray direction (or reflected if TIR occurs)
 */
Vector2 refract(const Vector2& incident, const Vector2& normal, float iorRatio, bool& totalInternalReflection) {
    float cosI = -dot(normal, incident);
    float sin2T = iorRatio * iorRatio * (1.0f - cosI * cosI);
    
    totalInternalReflection = false;
    
    if (sin2T >= 1.0f) {
        // Total internal reflection occurs
        totalInternalReflection = true;
        return reflect(incident, normal);
    }
    
    // Normal refraction using Snell's law
    float cosT = std::sqrt(1.0f - sin2T);
    return (iorRatio * incident) + (iorRatio * cosI - cosT) * normal;
}

/**
 * @brief Check if a point lies inside a triangle using barycentric coordinates
 * @param point Point to test
 * @param a First triangle vertex
 * @param b Second triangle vertex
 * @param c Third triangle vertex
 * @param epsilon Tolerance for edge cases
 * @return true if point is inside triangle, false otherwise
 */
bool pointInTriangle(const Vector2& point, const Vector2& a, const Vector2& b, const Vector2& c, float epsilon = 1e-6f) {
    Vector2 v0 = c - a;
    Vector2 v1 = b - a;
    Vector2 v2 = point - a;
    
    float dot00 = dot(v0, v0);
    float dot01 = dot(v0, v1);
    float dot02 = dot(v0, v2);
    float dot11 = dot(v1, v1);
    float dot12 = dot(v1, v2);
    
    // Compute barycentric coordinates
    float invDenom = 1.0f / (dot00 * dot11 - dot01 * dot01);
    float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    float v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    
    // Check if point is in triangle
    return (u >= -epsilon) && (v >= -epsilon) && (u + v <= 1.0f + epsilon);
}

/**
 * @brief Find where a ray first enters the prism
 * Handles intersection testing for all prism shapes with proper coordinate transformations
 * @param rayStart Ray starting point (world coordinates)
 * @param rayDir Ray direction vector (world coordinates)
 * @param cosTheta Cosine of prism rotation angle (pre-calculated for efficiency)
 * @param sinTheta Sine of prism rotation angle (pre-calculated for efficiency)
 * @param entryPoint Output: point where ray enters prism
 * @param entryNormal Output: surface normal at entry point
 * @param entryEdge Output: edge index where ray enters
 * @return true if ray hits prism, false otherwise
 */
bool findPrismEntry(const Vector2& rayStart, const Vector2& rayDir, 
                   float cosTheta, float sinTheta, // Pass pre-calculated values
                   Vector2& entryPoint, Vector2& entryNormal, int& entryEdge) {
    const float epsilon = 1e-4f;
    float bestT = 1e6f;
    bool found = false;
    
    // Transform ray to local prism coordinates (inverse rotation)
    
    Vector2 localRayStart(
        rayStart.x * cosTheta + rayStart.y * sinTheta,
        -rayStart.x * sinTheta + rayStart.y * cosTheta
    );
    
    Vector2 localRayDir(
        rayDir.x * cosTheta + rayDir.y * sinTheta,
        -rayDir.x * sinTheta + rayDir.y * cosTheta
    );
    
    switch (currentPrismType) {
        case RECTANGLE: {
            // Check intersection with each side
            Vector2 sides[4][2] = {
                { Vector2(-0.3f, -prismRectHeight), Vector2(0.3f, -prismRectHeight) }, // Bottom
                { Vector2(0.3f, -prismRectHeight), Vector2(0.3f, prismRectHeight) },   // Right
                { Vector2(0.3f, prismRectHeight), Vector2(-0.3f, prismRectHeight) },   // Top
                { Vector2(-0.3f, prismRectHeight), Vector2(-0.3f, -prismRectHeight) }  // Left
            };
            
            for (int i = 0; i < 4; i++) {
                Vector2 intersection;
                float t;
                if (raySegmentIntersection(localRayStart, localRayDir, sides[i][0], sides[i][1], intersection, t)) {
                    if (t > epsilon && t < bestT) {
                        bestT = t;
                        entryPoint = intersection;
                        entryEdge = i;
                        found = true;
                    }
                }
            }
            
            if (found) {
                // Calculate normal based on edge
                switch (entryEdge) {
                    case 0: entryNormal = Vector2(0.0f, -1.0f); break; // Bottom
                    case 1: entryNormal = Vector2(1.0f, 0.0f); break;  // Right
                    case 2: entryNormal = Vector2(0.0f, 1.0f); break;  // Top
                    case 3: entryNormal = Vector2(-1.0f, 0.0f); break; // Left
                }
            }
            break;
        }
        
        case TRIANGLE: {
            Vector2 vertices[3] = {
                Vector2(0.0f, 0.3f),   // Top
                Vector2(0.3f, -0.3f),  // Bottom right
                Vector2(-0.3f, -0.3f)  // Bottom left
            };
            
            for (int i = 0; i < 3; i++) {
                Vector2 segStart = vertices[i];
                Vector2 segEnd = vertices[(i + 1) % 3];
                Vector2 intersection;
                float t;
                
                if (raySegmentIntersection(localRayStart, localRayDir, segStart, segEnd, intersection, t)) {
                    if (t > epsilon && t < bestT) {
                        bestT = t;
                        entryPoint = intersection;
                        entryEdge = i;
                        found = true;
                    }
                }
            }
            
            if (found) {
                // Calculate normal (perpendicular to edge, pointing outward)
                Vector2 edge = vertices[(entryEdge + 1) % 3] - vertices[entryEdge];
                entryNormal = Vector2(-edge.y, edge.x).normalized();
                
                // Ensure normal points outward
                Vector2 center(0, -0.1f); // Approximate center of triangle
                if (dot(entryNormal, entryPoint - center) < 0) {
                    entryNormal = entryNormal * -1.0f;
                }
            }
            break;
        }
        
        case CIRCLE: {
            // Solve quadratic equation for ray-circle intersection
            float a = dot(localRayDir, localRayDir);
            float b = 2.0f * dot(localRayStart, localRayDir);
            float c = dot(localRayStart, localRayStart) - 0.3f * 0.3f;
            float discriminant = b * b - 4.0f * a * c;
            
            if (discriminant >= 0.0f) {
                float sqrtD = std::sqrt(discriminant);
                float t1 = (-b - sqrtD) / (2.0f * a);
                float t2 = (-b + sqrtD) / (2.0f * a);
                
                // Use nearest positive intersection
                float t = (t1 > epsilon) ? t1 : ((t2 > epsilon) ? t2 : -1.0f);
                if (t > epsilon) {
                    entryPoint = localRayStart + localRayDir * t;
                    entryNormal = entryPoint.normalized();
                    found = true;
                }
            }
            break;
        }
        
        case SEMICIRCLE: {
            bool foundArc = false, foundBase = false;
            Vector2 arcIntersect, baseIntersect;
            float arcT = 1e6f, baseT = 1e6f;
            
            // Check intersection with semicircular arc
            float a = dot(localRayDir, localRayDir);
            float b = 2.0f * dot(localRayStart, localRayDir);
            float c = dot(localRayStart, localRayStart) - 0.3f * 0.3f;
            float discriminant = b * b - 4.0f * a * c;
            
            if (discriminant >= 0.0f) {
                float sqrtD = std::sqrt(discriminant);
                float t1 = (-b - sqrtD) / (2.0f * a);
                float t2 = (-b + sqrtD) / (2.0f * a);
                
                // Check if intersections are in the upper half
                if (t1 > epsilon) {
                    Vector2 p1 = localRayStart + localRayDir * t1;
                    if (p1.y >= -epsilon) {
                        arcIntersect = p1;
                        arcT = t1;
                        foundArc = true;
                    }
                }
                
                if (t2 > epsilon && t2 < arcT) {
                    Vector2 p2 = localRayStart + localRayDir * t2;
                    if (p2.y >= -epsilon) {
                        arcIntersect = p2;
                        arcT = t2;
                        foundArc = true;
                    }
                }
            }
            
            // Check intersection with base line (y = 0)
            if (std::abs(localRayDir.y) > epsilon) {
                float t = -localRayStart.y / localRayDir.y;
                if (t > epsilon) {
                    float x = localRayStart.x + localRayDir.x * t;
                    if (x >= -0.3f - epsilon && x <= 0.3f + epsilon) {
                        baseIntersect = Vector2(x, 0.0f);
                        baseT = t;
                        foundBase = true;
                    }
                }
            }
            
            // Use the closest intersection
            if (foundArc && foundBase) {
                if (arcT < baseT) {
                    entryPoint = arcIntersect;
                    entryNormal = entryPoint.normalized();
                    entryEdge = 0; // Arc
                } else {
                    entryPoint = baseIntersect;
                    entryNormal = Vector2(0.0f, -1.0f);
                    entryEdge = 1; // Base
                }
                found = true;
            } else if (foundArc) {
                entryPoint = arcIntersect;
                entryNormal = entryPoint.normalized();
                entryEdge = 0; // Arc
                found = true;
            } else if (foundBase) {
                entryPoint = baseIntersect;
                // Normal should point outwards from the prism.
                // If ray hits from outside, normal is (0, -1).
                // If ray hits from inside, normal is (0, 1).
                // This depends on context not available here, but for entry it's usually from outside.
                entryNormal = Vector2(0.0f, -1.0f);
                entryEdge = 1; // Base
                found = true;
            }
            break;
        }
    }
    
    // Transform back to world coordinates
    if (found) {
        Vector2 worldEntryPoint(
            entryPoint.x * cosTheta - entryPoint.y * sinTheta,
            entryPoint.x * sinTheta + entryPoint.y * cosTheta
        );
        
        Vector2 worldEntryNormal(
            entryNormal.x * cosTheta - entryNormal.y * sinTheta,
            entryNormal.x * sinTheta + entryNormal.y * cosTheta
        );
        
        entryPoint = worldEntryPoint;
        entryNormal = worldEntryNormal;
    }
    
    return found;
}

// Find the *next* intersection of a ray that is *inside* the prism
bool findNextPrismIntersection(const Vector2& rayStart, const Vector2& rayDir,
                               int previousEdge,
                               float cosTheta, float sinTheta, // Pass pre-calculated values
                               Vector2& exitPoint, Vector2& exitNormal, int& exitEdge) {
    const float epsilon = 1e-4f;
    float bestT = 1e6f;
    bool found = false;

    // Transform ray to local coordinates
    Vector2 localRayStart(rayStart.x * cosTheta + rayStart.y * sinTheta, -rayStart.x * sinTheta + rayStart.y * cosTheta);
    Vector2 localRayDir(rayDir.x * cosTheta + rayDir.y * sinTheta, -rayDir.x * sinTheta + rayDir.y * cosTheta);

    // The logic is almost identical to findPrismEntry, but we must ignore the previous edge
    // and only accept intersections with t > epsilon.

    switch (currentPrismType) {
        case RECTANGLE: {
            Vector2 sides[4][2] = {
                { Vector2(-0.3f, -prismRectHeight), Vector2(0.3f, -prismRectHeight) },
                { Vector2(0.3f, -prismRectHeight), Vector2(0.3f, prismRectHeight) },
                { Vector2(0.3f, prismRectHeight), Vector2(-0.3f, prismRectHeight) },
                { Vector2(-0.3f, prismRectHeight), Vector2(-0.3f, -prismRectHeight) }
            };
            for (int i = 0; i < 4; i++) {
                if (i == previousEdge) continue; // Skip the edge we just came from
                Vector2 intersection;
                float t;
                if (raySegmentIntersection(localRayStart, localRayDir, sides[i][0], sides[i][1], intersection, t)) {
                    if (t > epsilon && t < bestT) {
                        bestT = t;
                        exitPoint = intersection;
                        exitEdge = i;
                        found = true;
                    }
                }
            }
            if (found) {
                switch (exitEdge) {
                    case 0: exitNormal = Vector2(0.0f, -1.0f); break;
                    case 1: exitNormal = Vector2(1.0f, 0.0f); break;
                    case 2: exitNormal = Vector2(0.0f, 1.0f); break;
                    case 3: exitNormal = Vector2(-1.0f, 0.0f); break;
                }
            }
            break;
        }
        case TRIANGLE: {
             Vector2 vertices[3] = { Vector2(0.0f, 0.3f), Vector2(0.3f, -0.3f), Vector2(-0.3f, -0.3f) };
             for (int i = 0; i < 3; i++) {
                if (i == previousEdge) continue;
                Vector2 segStart = vertices[i];
                Vector2 segEnd = vertices[(i + 1) % 3];
                Vector2 intersection;
                float t;
                if (raySegmentIntersection(localRayStart, localRayDir, segStart, segEnd, intersection, t)) {
                    if (t > epsilon && t < bestT) {
                        bestT = t;
                        exitPoint = intersection;
                        exitEdge = i;
                        found = true;
                    }
                }
            }
            if (found) {
                Vector2 edge = vertices[(exitEdge + 1) % 3] - vertices[exitEdge];
                exitNormal = Vector2(-edge.y, edge.x).normalized();
                Vector2 center(0, -0.1f);
                if (dot(exitNormal, exitPoint - center) < 0) {
                    exitNormal = exitNormal * -1.0f;
                }
            }
            break;
        }
        case CIRCLE: {
            float a = dot(localRayDir, localRayDir);
            float b = 2.0f * dot(localRayStart, localRayDir);
            float c = dot(localRayStart, localRayStart) - 0.3f * 0.3f;
            float discriminant = b * b - 4.0f * a * c;
            if (discriminant >= 0.0f) {
                float sqrtD = std::sqrt(discriminant);
                float t1 = (-b - sqrtD) / (2.0f * a);
                float t2 = (-b + sqrtD) / (2.0f * a);
                // We are inside, so we are looking for the next positive t value.
                // t1 is usually the point behind us, t2 is the one in front.
                if (t2 > epsilon) {
                    exitPoint = localRayStart + localRayDir * t2;
                    exitNormal = exitPoint.normalized();
                    found = true;
                }
            }
            break;
        }
        case SEMICIRCLE: {
            bool foundArc = false, foundBase = false;
            Vector2 arcIntersect, baseIntersect;
            float arcT = 1e6f, baseT = 1e6f;

            // Check arc (edge 0)
            if (previousEdge != 0) {
                float a = dot(localRayDir, localRayDir);
                float b = 2.0f * dot(localRayStart, localRayDir);
                float c = dot(localRayStart, localRayStart) - 0.3f * 0.3f;
                float discriminant = b * b - 4.0f * a * c;
                if (discriminant >= 0.0f) {
                    float t = (-b + std::sqrt(discriminant)) / (2.0f * a); // We are inside, so we need the positive t
                    if (t > epsilon) {
                        Vector2 p = localRayStart + localRayDir * t;
                        if (p.y >= -epsilon) {
                            arcIntersect = p;
                            arcT = t;
                            foundArc = true;
                        }
                    }
                }
            }

            // Check base (edge 1)
            if (previousEdge != 1) {
                if (std::abs(localRayDir.y) > epsilon) {
                    float t = -localRayStart.y / localRayDir.y;
                    if (t > epsilon) {
                        float x = localRayStart.x + localRayDir.x * t;
                        if (x >= -0.3f - epsilon && x <= 0.3f + epsilon) {
                            baseIntersect = Vector2(x, 0.0f);
                            baseT = t;
                            foundBase = true;
                        }
                    }
                }
            }

            if (foundArc && (!foundBase || arcT < baseT)) {
                exitPoint = arcIntersect;
                exitNormal = exitPoint.normalized();
                exitEdge = 0;
                found = true;
            } else if (foundBase) {
                exitPoint = baseIntersect;
                exitNormal = Vector2(0.0f, 1.0f); // Normal points outwards (up)
                exitEdge = 1;
                found = true;
            }
            break;
        }
    }

    if (found) {
        // Transform back to world coordinates
        Vector2 worldExitPoint(exitPoint.x * cosTheta - exitPoint.y * sinTheta, exitPoint.x * sinTheta + exitPoint.y * cosTheta);
        Vector2 worldExitNormal(exitNormal.x * cosTheta - exitNormal.y * sinTheta, exitNormal.x * sinTheta + exitNormal.y * cosTheta);
        exitPoint = worldExitPoint;
        exitNormal = worldExitNormal;
    }
    return found;
}


// Main simulation function
void draw_prism() {
    allRayAngles.clear();
    draw_prism_shape();

    // Pre-calculate rotation values
    float theta = prismRotation * M_PI / 180.0f;
    float cosTheta = std::cos(theta);
    float sinTheta = std::sin(theta);

    // 1. Initial Ray Setup
    float rayRadians = rotationAngle * M_PI / 180.0f;
    Vector2 initialRayDir(-std::cos(rayRadians), -std::sin(rayRadians));
    Vector2 initialRayStart = Vector2(0, 0) - initialRayDir * 2.0f; // Start far away

    // 2. Find First Entry Point
    Vector2 entryPoint, entryNormal;
    int entryEdge;
    if (!findPrismEntry(initialRayStart, initialRayDir, cosTheta, sinTheta, entryPoint, entryNormal, entryEdge)) {
        // Ray misses the prism, draw it to infinity
        Vector2 endPoint = initialRayStart + initialRayDir * 4.0f;
        draw_line(initialRayStart.x, initialRayStart.y, endPoint.x, endPoint.y,
                  rayAColor[0], rayAColor[1], rayAColor[2], 3.0f);
        return;
    }

    // Draw incident ray up to the prism
    draw_line(initialRayStart.x, initialRayStart.y, entryPoint.x, entryPoint.y,
              rayAColor[0], rayAColor[1], rayAColor[2], 3.0f);
    draw_point(entryPoint.x, entryPoint.y, 1, 1, 1, 10);

    // Calculate and display incident angle
    float cosI_entry = -dot(initialRayDir, entryNormal);
    rayIncidentAngle = std::acos(clamp(cosI_entry, -1, 1)) * 180.0f / M_PI;
    allRayAngles.push_back("Incident: " + std::to_string(rayIncidentAngle) + "°");

    // 3. First Refraction
    bool tir;
    Vector2 currentRayDir = refract(initialRayDir, entryNormal, refractiveIndex1 / refractiveIndex2, tir);
    currentRayDir.normalize();
    Vector2 currentRayStart = entryPoint + currentRayDir * 1e-4f; // Move slightly forward

    // 4. Trace Ray Inside Prism (with multiple bounces)
    const int MAX_BOUNCES = 10;
    for (int bounce = 0; bounce < MAX_BOUNCES; ++bounce) {
        // Find where the internal ray hits the prism boundary
        Vector2 exitPoint, exitNormal;
        int nextEdge; // Use a different name to avoid confusion
        if (!findNextPrismIntersection(currentRayStart, currentRayDir, entryEdge, cosTheta, sinTheta, exitPoint, exitNormal, nextEdge)) {
            // Ray is trapped or something went wrong. Draw to infinity as a fallback.
            Vector2 endPoint = currentRayStart + currentRayDir * 4.0f;
            draw_line(currentRayStart.x, currentRayStart.y, endPoint.x, endPoint.y,
                      rayBColor[0], rayBColor[1], rayBColor[2], 2.0f);
            break; // Exit loop
        }
        entryEdge = nextEdge; // Update the edge for the next iteration

        // Draw the internal ray segment
        draw_line(currentRayStart.x, currentRayStart.y, exitPoint.x, exitPoint.y,
                  rayBColor[0], rayBColor[1], rayBColor[2], 2.0f);
        draw_point(exitPoint.x, exitPoint.y, 1, 1, 1, 10);

        // 5. Decide to reflect (TIR) or refract at the exit point
        float n_in = refractiveIndex2;
        float n_out = refractiveIndex1;

        // The normal from findPrismEntry points outwards, which is what we need.
        float cosI = dot(currentRayDir, exitNormal);

        // If the ray and normal are pointing in the same general direction, something is wrong with the normal.
        // The normal should always oppose the incident ray direction for a correct cosI.
        // This indicates the normal is pointing inwards instead of outwards.
        if (cosI > 0) {
            // This case should ideally not be hit if normals are calculated correctly,
            // but as a safeguard, we flip the normal and recalculate cosI.
            exitNormal = exitNormal * -1.0f;
            cosI = dot(currentRayDir, exitNormal);
        }

        cosI = -cosI; // cosI should be positive

        float sin2T = (n_in / n_out) * (n_in / n_out) * (1.0f - cosI * cosI);

        if (sin2T >= 1.0f) { // Total Internal Reflection
            allRayAngles.push_back("TIR");
            currentRayDir = reflect(currentRayDir, exitNormal);
            currentRayDir.normalize();
            currentRayStart = exitPoint + currentRayDir * 1e-4f; // Move slightly forward
            // Continue to the next bounce in the loop
        } else { // Refract out of the prism
            float exitAngle = std::acos(clamp(cosI, -1, 1)) * 180.0f / M_PI;
            allRayAngles.push_back("Exit: " + std::to_string(exitAngle) + "°");

            Vector2 finalRayDir = refract(currentRayDir, exitNormal, n_in / n_out, tir);
            finalRayDir.normalize();
            Vector2 finalEndPoint = exitPoint + finalRayDir * 4.0f;
            draw_line(exitPoint.x, exitPoint.y, finalEndPoint.x, finalEndPoint.y,
                      rayCColor[0], rayCColor[1], rayCColor[2], 2.0f);
            break; // Ray has exited, so stop bouncing
        }
    }
}