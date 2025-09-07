/**
 * @file main.cpp
 * @brief Main application file for OpenGL Light Physics Simulation
 * 
 * This interactive educational tool demonstrates fundamental light physics including:
 * - Mirror reflection (law of reflection)
 * - Light refraction between media (Snell's law) 
 * - Total internal reflection
 * - Light dispersion and behavior in prisms
 * 
 * The application uses OpenGL for rendering, GLFW for window management,
 * GLEW for OpenGL extension loading, and ImGui for the user interface.
 */

#include "simulation.h"
#include <cmath>
#include <vector>

/**
 * @brief Supported prism shapes for the prism simulation mode
 */
enum PrismType {
    RECTANGLE,   ///< Rectangular prism
    TRIANGLE,    ///< Triangular prism  
    CIRCLE,      ///< Circular prism
    SEMICIRCLE   ///< Semicircular prism
};

// ============================================================================
// Global Variables Definitions
// ============================================================================

// Window and rendering properties
int screenWidth = 1080;                        ///< Initial window width
int screenHeight = 720;                        ///< Initial window height

// Simulation state
SimulationMode currentMode = REFLECTION;       ///< Current simulation mode
bool totalInternalReflection = false;         ///< Flag for total internal reflection
bool rayFromTop = false;                       ///< Direction flag for ray origin
float rotationAngle = 135.0f;                 ///< Current ray angle in degrees
float rayPosition = 0.0f;                     ///< Horizontal position of ray on surface

// Physics parameters
float refractiveIndex1 = 1.0f;                ///< First medium refractive index (air)
float refractiveIndex2 = 1.33f;               ///< Second medium refractive index (water)
float n1 = 1.0f, n2 = 1.33f;                  ///< Active refractive indices for calculations

// Visual customization
float rayAColor[3] = { 1.0f, 0.8f, 0.0f };    ///< Incident ray color (golden yellow)
float rayBColor[3] = { 0.0f, 0.8f, 1.0f };    ///< Refracted ray color (cyan blue)
float rayCColor[3] = { 0.4f, 0.1f, 0.5f };    ///< Reflected ray color (purple)
float refractiveIndex2_RGB[3] = { 1.0f, 1.0f, 1.0f }; ///< RGB-specific refractive indices

// Mouse interaction state
double lastX = 0.0;                           ///< Last mouse X position for tracking
const float interfaceY = 0.0f;               ///< Y coordinate of refraction interface
bool isFollowingMouse = false;                ///< Flag for mouse angle control mode
bool isDraggingRay = false;                   ///< Flag for ray position dragging mode

// ============================================================================
// Callback and Utility Functions
// ============================================================================

/**
 * @brief Framebuffer size callback for window resizing
 * Updates the OpenGL viewport and projection matrix when window is resized
 * @param window GLFW window handle
 * @param width New window width
 * @param height New window height
 */
void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    screenWidth = width;
    screenHeight = height;
    glViewport(0, 0, width, height);
    
    // Set up orthographic projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float aspect = (float)width / (float)height;
    
    // Maintain aspect ratio while covering the full screen
    if (aspect >= 1.0f) {
        glOrtho(-aspect, aspect, -1.0f, 1.0f, -1.0f, 1.0f);
    } else {
        glOrtho(-1.0f, 1.0f, -1.0f / aspect, 1.0f / aspect, -1.0f, 1.0f);
    }
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

/**
 * @brief Determine if ray originates from top based on angle
 * @param angle Ray angle in degrees
 */
void isRayFromTop(float angle) {
    rayFromTop = (angle <= 180.0f);
}

/**
 * @brief Mouse button callback for interactive ray control
 * Handles two interaction modes:
 * 1. Click and drag to change ray angle (mouse following mode)
 * 2. Click and drag ray intersection point to move ray position
 * @param window GLFW window handle
 * @param button Mouse button (left, right, middle)
 * @param action Press or release action
 * @param mods Modifier keys (shift, ctrl, alt)
 */
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    ImGuiIO& io = ImGui::GetIO();
    // Don't process mouse input if ImGui wants to capture it
    if (io.WantCaptureMouse)
        return;
    
    // Get current mouse position and convert to scene coordinates
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);
    float normX = static_cast<float>(xpos) / screenWidth * 2.0f - 1.0f;
    float normY = 1.0f - static_cast<float>(ypos) / screenHeight * 2.0f;
    float aspect = (float)screenWidth / (float)screenHeight;
    float sceneX, sceneY;
    
    // Apply aspect ratio correction
    if (aspect >= 1.0f) {
        sceneX = normX * aspect;
        sceneY = normY;
    } else {
        sceneX = normX;
        sceneY = normY / aspect;
    }
    
    // Get surface Y coordinate based on current simulation mode
    float surfaceY = (currentMode == REFLECTION) ? -0.8f : interfaceY;
    
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            // Check if clicking near the ray intersection point for dragging
            if (std::fabs(sceneX - rayPosition) < 0.05f && std::fabs(sceneY - surfaceY) < 0.05f) {
                isDraggingRay = true;
            } else {
                // Start mouse following mode for angle control
                isFollowingMouse = true;
                float dx = sceneX - rayPosition;
                float dy = sceneY - surfaceY;
                float angleRad = std::atan2(dy, dx);
                float angleDeg = angleRad * 180.0f / M_PI;
                if(angleDeg < 0) angleDeg += 360.0f;
                rotationAngle = angleDeg;
                
                // Constrain angle for reflection mode
                if (currentMode == REFLECTION) {
                    if (rotationAngle < 0.0f) rotationAngle = 0.0f;
                    if (rotationAngle > 180.0f) rotationAngle = 180.0f;
                }
            }
        }
        else if (action == GLFW_RELEASE) {
            // End both interaction modes
            isFollowingMouse = false;
            isDraggingRay = false;
        }
    }
}

/**
 * @brief Cursor position callback for continuous mouse interaction
 * Updates ray angle or position based on current interaction mode
 * @param window GLFW window handle
 * @param xpos Mouse X position in screen coordinates
 * @param ypos Mouse Y position in screen coordinates
 */
void cursor_position_callback(GLFWwindow* window, double xpos, double ypos) {
    ImGuiIO& io = ImGui::GetIO();
    // Don't process mouse input if ImGui wants to capture it
    if (io.WantCaptureMouse)
        return;
    
    // Convert screen coordinates to normalized device coordinates
    float normX = static_cast<float>(xpos) / screenWidth * 2.0f - 1.0f;
    float normY = 1.0f - static_cast<float>(ypos) / screenHeight * 2.0f;
    float aspect = (float)screenWidth / (float)screenHeight;
    float sceneX, sceneY;
    
    // Apply aspect ratio correction to get scene coordinates
    if (aspect >= 1.0f) {
        sceneX = normX * aspect;
        sceneY = normY;
    } else {
        sceneX = normX;
        sceneY = normY / aspect;
    }
    
    float surfaceY = (currentMode == REFLECTION) ? -0.8f : interfaceY;
    
    if (isFollowingMouse) {
        // Update ray angle based on mouse position relative to ray intersection
        float dx = sceneX - rayPosition;
        float dy = sceneY - surfaceY;
        float angleRad = std::atan2(dy, dx);
        float angleDeg = angleRad * 180.0f / M_PI;
        rotationAngle = angleDeg;
        
        // Normalize angle to 0-360 range
        rotationAngle = fmod(rotationAngle, 360.0f);
        if (rotationAngle < 0) rotationAngle += 360.0f;
        
        // Apply constraints for reflection mode
        if (currentMode == REFLECTION) {
            if (rotationAngle < 0.0f) rotationAngle = 0.0f;
            if (rotationAngle > 180.0f) rotationAngle = 180.0f;
        }
    }
    else if (isDraggingRay) {
        // Update ray position by tracking horizontal mouse movement
        double deltaX = xpos - lastX;
        rayPosition += static_cast<float>(deltaX) * 0.005f;
        
        // Clamp ray position to visible area
        if (rayPosition < -0.9f) rayPosition = -0.9f;
        if (rayPosition > 0.9f) rayPosition = 0.9f;
    }
    lastX = xpos;
}

// ============================================================================
// Main Application Function
// ============================================================================

/**
 * @brief Main function - application entry point
 * Initializes OpenGL, GLFW, GLEW, and ImGui, then runs the main render loop
 * @return 0 on success, -1 on failure
 */
int main(){
    // Initialize GLFW
    if (!glfwInit())
        return -1;
    
    // Create window
    GLFWwindow* window = glfwCreateWindow(screenWidth, screenHeight, "Light Simulation", NULL, NULL);
    if (!window) {
        glfwTerminate();
        return -1;
    }
    
    // Set up OpenGL context
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);  // Enable v-sync
    
    // Initialize GLEW for OpenGL extension loading
    if (glewInit() != GLEW_OK) {
        std::cerr << "Error initializing GLEW" << std::endl;
        glfwTerminate();
        return -1;
    }
    
    // Set up GLFW callbacks for user interaction
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    framebuffer_size_callback(window, screenWidth, screenHeight);
    
    // Initialize ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");
    
    // ========================================================================
    // Main Render Loop
    // ========================================================================
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        
        // Start new ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();
        
        // ====================================================================
        // User Interface (ImGui Controls Window)
        // ====================================================================
        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(300, 450), ImGuiCond_FirstUseEver);
        ImGui::Begin("Controls");
        
        // Simulation mode selection
        const char* modes[] = { "Reflection", "Refraction", "Prism" };
        static int currentModeIndex = static_cast<int>(currentMode);
        if (ImGui::Combo("Simulation Mode", &currentModeIndex, modes, IM_ARRAYSIZE(modes))) {
            currentMode = static_cast<SimulationMode>(currentModeIndex);
        }
        
        ImGui::Separator();
        
        // Mode-specific controls
        if (currentMode == REFLECTION) {
            // Reflection mode controls
            ImGui::Text("Ray Angle");
            ImGui::SliderFloat("##AngleSlider", &rotationAngle, 0.0f, 180.0f, "%.1f°");
            float incident = std::fabs(90.0f - (180.0f - rotationAngle));
            ImGui::Separator();
            ImGui::Text("Incident Angle: %.1f°", incident);
            ImGui::Text("Reflection Angle: %.1f°", incident);
        }
        else if (currentMode == REFRACTION) {
            // Refraction mode controls
            isRayFromTop(rotationAngle);
            ImGui::Text("Ray Angle");
            ImGui::SliderFloat("##AngleSlider", &rotationAngle, 0.0f, 360.0f, "%.1f°");
            
            ImGui::Separator();
            ImGui::Text("Refractive Indices");
            ImGui::SliderFloat("Top Medium (n1)", &refractiveIndex1, 1.0f, 3.0f, "%.2f");
            ImGui::SliderFloat("Bottom Medium (n2)", &refractiveIndex2, 1.0f, 3.0f, "%.2f");
            
            // Update active refractive indices based on ray direction
            n1 = rayFromTop ? refractiveIndex1 : refractiveIndex2;
            n2 = rayFromTop ? refractiveIndex2 : refractiveIndex1;
            ImGui::Text("From: n = %.2f  To: n = %.2f", n1, n2);
            
            // Calculate and display angles
            float incident = std::fabs(90.0f - std::fmod(rotationAngle, 360.0f));
            ImGui::Text("Incident Angle: %.1f°", incident);
            
            // Check for total internal reflection using Snell's law
            float incidentRad = incident * M_PI / 180.0f;
            float sinRefr = (n1 / n2) * sinf(incidentRad);

            if (sinRefr > 1.0f || sinRefr < -1.0f) {
                totalInternalReflection = true;
            } else {
                float refrAngle = asinf(sinRefr) * 180.0f / M_PI;
                totalInternalReflection = false;
            }
            
            // Display refraction angle or total internal reflection warning
            if (!totalInternalReflection) {
                float refrAngle = asinf(sinRefr) * 180.0f / M_PI;
                ImGui::Text("Refraction Angle: %.1f°", refrAngle);
            } else {
                ImGui::TextColored(ImVec4(1.0f, 0.3f, 0.3f, 1.0f), "Total Internal Reflection");
            }
        }
        else if (currentMode == PRISM) {
            // Prism mode controls
            const char* prismas[] = { "Rectangle", "Triangle", "Circle", "Semicircle" }; 
            extern PrismType currentPrismType;
            static int currentPrismIndex = static_cast<int>(currentPrismType);
            ImGui::Combo("Prism shape", &currentPrismIndex, prismas, IM_ARRAYSIZE(prismas));
            if (static_cast<PrismType>(currentPrismIndex) != currentPrismType) {
                currentPrismType = static_cast<PrismType>(currentPrismIndex);
            }
            
            ImGui::SliderFloat("Incident ray", &rotationAngle, 0.0f, 360.0f, "%.1f°");
            extern float prismRotation;
            ImGui::SliderFloat("Prism Rotation", &prismRotation, 0.0f, 360.0f, "%.1f°");
            
            // Rectangle height slider (only for rectangular prisms)
            extern float prismRectHeight;
            if (currentPrismType == RECTANGLE) {
                ImGui::SliderFloat("Rectangle Height", &prismRectHeight, 0.05f, 0.98f, "%.2f");
            }
            
            ImGui::Separator();
            // Display ray angle information
            extern float rayIncidentAngle, rayInternalAngle, rayExitAngle, rayDistance;
            ImGui::Text("Incident Ray Angle:   %.2f°", rayIncidentAngle);
            ImGui::Text("Internal Ray Angle:   %.2f°", rayInternalAngle);
            ImGui::Text("Exit Ray Angle:       %.2f°", rayExitAngle);
            ImGui::Text("Ray Distance: %.3f", rayDistance);
            
            ImGui::Separator();
            ImGui::Text("Refractive Indices");
            ImGui::SliderFloat("Environmental Medium (n1)", &refractiveIndex1, 1.0f, 3.0f, "%.2f");
            ImGui::SliderFloat3("Prism Medium (n2) RGB", refractiveIndex2_RGB, 1.0f, 3.0f, "%.2f");
            
            // Display all ray angles for complex prism interactions
            extern std::vector<std::string> allRayAngles;
            ImGui::Text("Ray Angles (all events):");
            for (const auto& angleStr : allRayAngles) {
                ImGui::TextUnformatted(angleStr.c_str());
            }
            ImGui::Separator();
        }
        
        // Ray color customization
        ImGui::Separator();
        ImGui::Text("Ray Colors");
        ImGui::ColorEdit3("Incident Ray (A)", rayAColor);
        if (currentMode == REFLECTION) {
            ImGui::ColorEdit3("Reflected Ray (B)", rayBColor);
        } else if(currentMode == REFRACTION || currentMode == PRISM) {
            ImGui::ColorEdit3("Refracted Ray (B)", rayBColor);
            ImGui::ColorEdit3("Reflected Ray (C)", rayCColor);
        }
        ImGui::End();
        
        // ====================================================================
        // 3D Scene Rendering
        // ====================================================================
        glClear(GL_COLOR_BUFFER_BIT);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        
        // Draw simulation based on active mode
        if (currentMode == REFLECTION) {
            draw_reflection();
        } else if (currentMode == REFRACTION) {
            draw_refraction();
        } else if (currentMode == PRISM) {
            draw_prism();
        }
        
        // Render ImGui interface
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }
    
    // ========================================================================
    // Cleanup
    // ========================================================================
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}