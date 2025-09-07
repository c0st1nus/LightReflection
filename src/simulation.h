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
#include <cstdio>  
#include <iostream>  

// ============================================================================
// Global Variables (definitions in main.cpp)
// ============================================================================

// Screen dimensions
extern int screenWidth;  
extern int screenHeight;  

// Simulation modes
enum SimulationMode { REFLECTION, REFRACTION, PRISM };  
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

#endif // SIMULATION_H