/**
 * @file simulation.h
 * @brief Main header file for the OpenGL Light Physics Simulation
 * 
 * This file contains all the global declarations, enums, and utility functions
 * needed for simulating light behavior including reflection, refraction, and
 * prism interactions.
 */

#ifndef SIMULATION_H  
#define SIMULATION_H  

#include <GL/glew.h>  
#include <GLFW/glfw3.h>  
#define _USE_MATH_DEFINES
#include <cmath>  
#include "imgui.h"  
#include "imgui_impl_glfw.h"  
#include "imgui_impl_opengl3.h"  
#include <algorithm>
#include <cstdio>  
#include <iostream>  
#include <string>
#include "mirror.h"
#include "geometric_optics.h"

// ============================================================================
// Global Variables (definitions in main.cpp)
// ============================================================================

// Screen dimensions
extern int screenWidth;  
extern int screenHeight;  

// Simulation modes
enum SimulationMode { REFLECTION, REFRACTION, PRISM, MIRROR, LENS };  
extern SimulationMode currentMode;  

// Physics state variables
extern bool totalInternalReflection;    // Flag for total internal reflection
extern bool rayFromTop;                 // Direction flag for ray origin
extern float rotationAngle;             // Current ray angle in degrees
extern float rayPosition;               // Horizontal position of ray on surface

// Refractive indices for different media
extern float refractiveIndex1;          // First medium (usually air = 1.0)
extern float refractiveIndex2;          // Second medium (glass, water, etc.)
extern float n1, n2;                    // Alternative names for refractive indices
extern float refractiveIndex2_RGB[3];   // RGB-specific refractive indices for dispersion

// Ray color customization (RGB values 0.0-1.0)
extern float rayAColor[3];              // Incident ray color
extern float rayBColor[3];              // Transmitted/refracted ray color  
extern float rayCColor[3];              // Reflected ray color  

// Object-mode optical configuration
extern MirrorType currentMirrorType;
extern float mirrorCurvatureRadius;
extern LensType currentLensType;
extern float lensFocalLength;
extern float objectDistance;
extern float objectHeight;
extern bool objectPointsUp;
extern bool showConstructionRays;

// ============================================================================
// Utility Drawing Functions
// ============================================================================

/**
 * @brief Draw a point at specified coordinates
 * @param x X coordinate in normalized device coordinates
 * @param y Y coordinate in normalized device coordinates  
 * @param r Red component (0.0-1.0)
 * @param g Green component (0.0-1.0)
 * @param b Blue component (0.0-1.0)
 * @param size Point size in pixels
 */
inline void draw_point(float x, float y, float r, float g, float b, float size) {  
   GLfloat point[] = { x, y };  
   glEnableClientState(GL_VERTEX_ARRAY);  
   glPointSize(size);  
   glColor3f(r, g, b);  
   glVertexPointer(2, GL_FLOAT, 0, point);  
   glDrawArrays(GL_POINTS, 0, 1);  
   glDisableClientState(GL_VERTEX_ARRAY);  
}  

/**
 * @brief Draw a solid line between two points
 * @param x1 Start point X coordinate
 * @param y1 Start point Y coordinate
 * @param x2 End point X coordinate
 * @param y2 End point Y coordinate
 * @param r Red component (0.0-1.0)
 * @param g Green component (0.0-1.0)
 * @param b Blue component (0.0-1.0)
 * @param width Line width in pixels
 */
inline void draw_line(float x1, float y1, float x2, float y2, float r, float g, float b, float width) {  
   glLineWidth(width);  
   glBegin(GL_LINES);  
   glColor3f(r, g, b);  
   glVertex2f(x1, y1);  
   glVertex2f(x2, y2);  
   glEnd();  
}  

/**
 * @brief Draw a dashed line between two points
 * @param x1 Start point X coordinate
 * @param y1 Start point Y coordinate
 * @param x2 End point X coordinate
 * @param y2 End point Y coordinate
 * @param r Red component (0.0-1.0)
 * @param g Green component (0.0-1.0)
 * @param b Blue component (0.0-1.0)
 * @param width Line width in pixels
 */
inline void draw_dashed_line(float x1, float y1, float x2, float y2, float r, float g, float b, float width) {  
   glLineStipple(2, 0xAAAA);  
   glEnable(GL_LINE_STIPPLE);  
   glLineWidth(width);  
   glBegin(GL_LINES);  
   glColor3f(r, g, b);  
   glVertex2f(x1, y1);  
   glVertex2f(x2, y2);  
   glEnd();  
   glDisable(GL_LINE_STIPPLE);  
}  

inline void draw_filled_triangle(const Vector2& a, const Vector2& b, const Vector2& c, float r, float g, float bColor) {
   glBegin(GL_TRIANGLES);
   glColor3f(r, g, bColor);
   glVertex2f(a.x, a.y);
   glVertex2f(b.x, b.y);
   glVertex2f(c.x, c.y);
   glEnd();
}

inline float scene_half_width() {
   float aspect = (screenHeight > 0) ? static_cast<float>(screenWidth) / static_cast<float>(screenHeight) : 1.0f;
   return (aspect >= 1.0f) ? aspect : 1.0f;
}

inline float scene_half_height() {
   float aspect = (screenHeight > 0) ? static_cast<float>(screenWidth) / static_cast<float>(screenHeight) : 1.0f;
   return (aspect >= 1.0f) ? 1.0f : (1.0f / aspect);
}

inline Vector2 screen_to_scene(double xpos, double ypos) {
   float normX = static_cast<float>(xpos) / static_cast<float>(screenWidth) * 2.0f - 1.0f;
   float normY = 1.0f - static_cast<float>(ypos) / static_cast<float>(screenHeight) * 2.0f;
   return Vector2(normX * scene_half_width(), normY * scene_half_height());
}

inline ImVec2 scene_to_screen(const Vector2& point) {
   float halfWidth = scene_half_width();
   float halfHeight = scene_half_height();
   float screenX = ((point.x + halfWidth) / (2.0f * halfWidth)) * static_cast<float>(screenWidth);
   float screenY = ((halfHeight - point.y) / (2.0f * halfHeight)) * static_cast<float>(screenHeight);
   return ImVec2(screenX, screenY);
}

inline float diagram_axis_min_x() {
   float halfWidth = scene_half_width();
   return -halfWidth + 0.10f * halfWidth;
}

inline float diagram_axis_max_x() {
   float halfWidth = scene_half_width();
   return halfWidth - 0.10f * halfWidth;
}

inline float mirror_diagram_x() {
   float minX = diagram_axis_min_x();
   float maxX = diagram_axis_max_x();
   return minX + 0.80f * (maxX - minX);
}

inline float lens_diagram_x() {
   return 0.0f;
}

inline float diagram_object_height_limit() {
   return 0.55f * scene_half_height();
}

inline float mirror_object_distance_limit() {
   return mirror_diagram_x() - (diagram_axis_min_x() + 0.12f * scene_half_width());
}

inline float lens_object_distance_limit() {
   return lens_diagram_x() - (diagram_axis_min_x() + 0.12f * scene_half_width());
}

inline float clamp_mirror_object_distance(float distance) {
   return std::clamp(distance, 0.18f, mirror_object_distance_limit());
}

inline float clamp_lens_object_distance(float distance) {
   return std::clamp(distance, 0.18f, lens_object_distance_limit());
}

inline void draw_scene_label(const Vector2& point, const std::string& text, const ImVec4& color, bool centered = true) {
   ImDrawList* drawList = ImGui::GetBackgroundDrawList();
   ImVec2 screen = scene_to_screen(point);
   ImVec2 textSize = ImGui::CalcTextSize(text.c_str());
   if (centered) {
      screen.x -= textSize.x * 0.5f;
   }
   drawList->AddText(screen, ImGui::ColorConvertFloat4ToU32(color), text.c_str());
}

// ============================================================================
// Function Prototypes for Simulation-Specific Drawing Functions
// ============================================================================

/**
 * @brief Draw the reflection simulation scene
 * Shows incident and reflected rays bouncing off a mirror surface
 */
void draw_reflection();  

/**
 * @brief Draw the refraction simulation scene  
 * Shows light bending when passing between different media
 */
void draw_refraction();  

/**
 * @brief Draw the prism simulation scene
 * Shows light dispersion and internal reflections in geometric prisms
 */
void draw_prism();  

/**
 * @brief Draw the object-based mirror diagram scene
 */
void draw_mirror_mode();

/**
 * @brief Draw the object-based lens diagram scene
 */
void draw_lens_mode();

#endif // SIMULATION_H
