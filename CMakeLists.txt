# ============================================================================
# OpenGL Light Physics Simulation - CMake Build Configuration with vcpkg
# ============================================================================

cmake_minimum_required(VERSION 3.15)

# Bootstrap a local vcpkg checkout when no toolchain was provided explicitly.
if(NOT DEFINED CMAKE_TOOLCHAIN_FILE)
    set(VCPKG_ROOT "${CMAKE_SOURCE_DIR}/.vcpkg" CACHE PATH "Local vcpkg checkout")
    set(_vcpkg_toolchain "${VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake")

    if(NOT EXISTS "${_vcpkg_toolchain}")
        find_package(Git REQUIRED)

        message(STATUS "vcpkg toolchain not provided. Bootstrapping local vcpkg into ${VCPKG_ROOT}")

        execute_process(
            COMMAND "${GIT_EXECUTABLE}" clone https://github.com/microsoft/vcpkg.git "${VCPKG_ROOT}"
            COMMAND_ERROR_IS_FATAL ANY
        )

        if(WIN32)
            execute_process(
                COMMAND "${VCPKG_ROOT}/bootstrap-vcpkg.bat" -disableMetrics
                WORKING_DIRECTORY "${VCPKG_ROOT}"
                COMMAND_ERROR_IS_FATAL ANY
            )
        else()
            execute_process(
                COMMAND "${CMAKE_COMMAND}" -E env VCPKG_DISABLE_METRICS=1 sh "${VCPKG_ROOT}/bootstrap-vcpkg.sh" -disableMetrics
                WORKING_DIRECTORY "${VCPKG_ROOT}"
                COMMAND_ERROR_IS_FATAL ANY
            )
        endif()
    endif()

    set(CMAKE_TOOLCHAIN_FILE "${_vcpkg_toolchain}" CACHE STRING "vcpkg toolchain file")
endif()

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
    src/geometric_optics.cpp # Thin-lens / mirror image calculations
    src/main.cpp           # Main application and UI
    src/mirror.cpp         # Mirror geometry and curved mirror support
    src/mirror_mode.cpp    # Object-based mirror diagram mode
    src/lens_mode.cpp      # Object-based lens diagram mode
    src/optics.cpp         # Shared optical geometry and tracing helpers
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
