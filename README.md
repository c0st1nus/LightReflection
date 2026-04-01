# OpenGL Light Physics Simulation

An interactive educational tool that demonstrates fundamental principles of light physics through real-time OpenGL visualization.

![Light Simulation Demo](docs/demo.gif)

## Overview

This application provides three interactive simulation modes to explore different aspects of light behavior:

1. **Reflection** - Mirror reflection following the law of reflection
2. **Refraction** - Light bending between different media using Snell's law
3. **Prism** - Complex light interactions with geometric prisms

## Features

### 🔍 **Reflection Mode**
- Interactive ray angle control
- Real-time calculation of incident and reflection angles
- Visual demonstration that angle of incidence equals angle of reflection
- Mouse-controlled ray positioning and direction

### 🌊 **Refraction Mode**
- Adjustable refractive indices for different media (air, water, glass, etc.)
- Snell's law calculations with real-time angle display
- Total internal reflection visualization
- Bidirectional ray simulation (top-to-bottom and bottom-to-top)

### 💎 **Prism Mode**
- Multiple prism shapes: Rectangle, Triangle, Circle, Semicircle
- Complex ray tracing with multiple internal reflections
- RGB dispersion simulation with wavelength-dependent refractive indices
- Fresnel reflection coefficients
- Ray path tracking and angle logging

### 🎨 **Interactive Controls**
- Real-time ray color customization
- Mouse-driven ray angle and position control
- Adjustable prism rotation and geometry
- Live physics parameter feedback

## Technical Details

### Built With
- **OpenGL** - 2D graphics rendering
- **GLFW** - Window management and input handling  
- **GLEW** - OpenGL extension loading
- **ImGui** - Immediate mode GUI for controls
- **C++17** - Modern C++ with mathematical computations

### Physics Algorithms
- **Snell's Law**: `n₁ sin(θ₁) = n₂ sin(θ₂)` for refraction calculations
- **Fresnel Equations**: For reflection/transmission coefficients
- **Vector Mathematics**: 2D ray-surface intersection algorithms
- **Geometric Computations**: Point-in-polygon testing for complex prisms

## Building and Installation

### Prerequisites
- **C++17 compatible compiler** (MSVC 2019+, GCC 8+, or Clang 9+)
- **CMake 3.15+**
- **Git** (used by CMake to bootstrap `vcpkg`)

### Dependency Management
The project uses `vcpkg` manifest mode via [`vcpkg.json`](vcpkg.json).

On the first CMake configure, [`CMakeLists.txt`](CMakeLists.txt) will:
- clone `vcpkg` into `.vcpkg/` if no toolchain was provided manually
- bootstrap the local `vcpkg` executable
- install the required dependencies declared in `vcpkg.json`

The dependency baseline is pinned, so fresh checkouts resolve the same package set by default.

### Build Instructions

#### Windows (Visual Studio)
```bash
git clone https://github.com/yourusername/OpenGL-Cmake.git
cd OpenGL-Cmake
cmake -S . -B build -G "Visual Studio 16 2019" -A x64
cmake --build build --config Release
```

#### Windows (MinGW)
```bash
git clone https://github.com/yourusername/OpenGL-Cmake.git
cd OpenGL-Cmake
cmake -S . -B build -G "MinGW Makefiles"
cmake --build build
```

#### Linux
```bash
git clone https://github.com/yourusername/OpenGL-Cmake.git
cd OpenGL-Cmake
cmake -S . -B build
cmake --build build -j$(nproc)
```

#### macOS
```bash
git clone https://github.com/yourusername/OpenGL-Cmake.git
cd OpenGL-Cmake
cmake -S . -B build
cmake --build build -j$(sysctl -n hw.ncpu)
```

#### Optional: Reuse an Existing vcpkg Installation
If you already have `vcpkg`, you can keep using it instead of the auto-bootstrapped `.vcpkg/` checkout:

```bash
cmake -S . -B build \
  -DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake
cmake --build build
```

### Running the Application
```bash
# Windows
.\build\Release\opengl_simulation.exe

# Linux/macOS
./build/opengl_simulation
```

## Usage Guide

### Basic Controls
- **Left Click + Drag**: Change ray angle by pointing mouse cursor
- **Left Click + Drag** (near intersection): Move ray position along surface
- **Mode Selection**: Choose between Reflection, Refraction, or Prism modes
- **Parameter Sliders**: Adjust refractive indices, angles, and colors in real-time

### Tips for Exploration
1. **Start with Reflection**: Observe how changing the ray angle affects reflection
2. **Try Total Internal Reflection**: In refraction mode, set n₁ > n₂ and increase the angle
3. **Experiment with Prisms**: Rotate different prism shapes and observe light paths
4. **Custom Materials**: Use realistic refractive indices:
   - Air: 1.00
   - Water: 1.33
   - Glass: 1.5-1.9
   - Diamond: 2.42

## Project Structure

```
OpenGL-Cmake/
├── CMakeLists.txt          # Build configuration
├── src/
│   ├── main.cpp            # Application entry point and UI
│   ├── simulation.h        # Global declarations and utility functions
│   ├── vector2.h           # 2D vector mathematics
│   ├── ray.h               # Ray class header
│   ├── ray.cpp             # Ray class implementation
│   ├── reflection.cpp      # Mirror reflection simulation
│   ├── refraction.cpp      # Light refraction simulation
│   ├── prism.cpp           # Complex prism light interactions
├── glfw/                  # GLFW library (Windows binaries)
├── glew/                  # GLEW library (Windows binaries)
├── imgui/                 # ImGui source files
└── build-*/               # Build directories (generated)
```

## Educational Applications

This simulation is ideal for:
- **Physics Education**: Visualizing abstract optical concepts
- **Interactive Learning**: Hands-on exploration of light behavior
- **STEM Demonstrations**: Real-time physics parameter manipulation
- **Computer Graphics**: Understanding ray tracing fundamentals

## Contributing

Contributions are welcome! Areas for improvement:
- Additional prism shapes (hexagon, pentagon, custom polygons)
- Chromatic dispersion visualization
- Wave optics phenomena (interference, diffraction)
- Performance optimizations for complex ray tracing
- Educational mode with guided tutorials

### Development Setup
1. Fork the repository
2. Create a feature branch: `git checkout -b feature/new-feature`
3. Make changes and test thoroughly
4. Commit with descriptive messages
5. Push and create a pull request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **GLFW Team** - Cross-platform window management
- **GLEW Team** - OpenGL extension loading
- **ImGui Team** - Immediate mode GUI framework
- **OpenGL Community** - Graphics programming resources
- **Physics References** - Optics textbooks and online resources

## Physics References

- Hecht, Eugene. "Optics" (5th Edition)
- Born, Max & Wolf, Emil. "Principles of Optics"
- [Snell's Law on Wikipedia](https://en.wikipedia.org/wiki/Snell%27s_law)
- [Fresnel Equations on Wikipedia](https://en.wikipedia.org/wiki/Fresnel_equations)

---

**Note**: This is an educational simulation. While physically accurate for geometric optics, it simplifies some aspects of real light behavior for clarity and performance.
