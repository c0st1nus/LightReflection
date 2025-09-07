# ============================================================================
# OpenGL Light Physics Simulation - CMake Build Configuration with vcpkg
# ============================================================================

cmake_minimum_required(VERSION 3.15)
project(OpenGLSimulation LANGUAGES CXX)

# ============================================================================
# Compiler Configuration
# ============================================================================

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# ============================================================================
# Find packages installed by vcpkg
# ============================================================================

find_package(GLEW REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(imgui CONFIG REQUIRED)
find_package(OpenGL REQUIRED)

# ============================================================================
# Source Files Configuration
# ============================================================================

# Application source files
set(APP_SOURCES
    src/main.cpp           # Main application and UI
    src/prism.cpp          # Complex prism light simulation
    src/ray.cpp            # Ray class implementation
    src/reflection.cpp     # Mirror reflection simulation
    src/refraction.cpp     # Light refraction simulation
)

# ============================================================================
# Executable Configuration
# ============================================================================

# Create executable target
add_executable(opengl_simulation ${APP_SOURCES})

# ============================================================================
# Link Libraries (vcpkg packages)
# ============================================================================

target_link_libraries(opengl_simulation PRIVATE 
    GLEW::GLEW              # GLEW OpenGL extension loading
    glfw                    # GLFW window management
    imgui::imgui            # ImGui immediate mode GUI
    OpenGL::GL              # OpenGL library
)

# For Windows, may need additional system libraries
if(WIN32)
    target_link_libraries(opengl_simulation PRIVATE 
        gdi32               # Windows Graphics Device Interface
    )
endif()
