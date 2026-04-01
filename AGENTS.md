# Repository Guidelines

## Project Structure & Module Organization
Core application code lives in [`src/`](./src). [`src/main.cpp`](./src/main.cpp) owns window setup, ImGui integration, and mode switching. Physics behavior is split by feature: [`reflection.cpp`](./src/reflection.cpp), [`refraction.cpp`](./src/refraction.cpp), and [`prism.cpp`](./src/prism.cpp). Shared math and ray primitives live in [`vector2.h`](./src/vector2.h), [`ray.h`](./src/ray.h), and [`ray.cpp`](./src/ray.cpp). Build configuration is in [`CMakeLists.txt`](./CMakeLists.txt); dependency pinning is in [`vcpkg.json`](./vcpkg.json).

## Build, Test, and Development Commands
- `cmake -S . -B build`: configure the project. If no toolchain is provided, CMake bootstraps local `vcpkg` into `.vcpkg/` and installs dependencies automatically.
- `cmake --build build -j$(nproc)`: compile `opengl_simulation`.
- `./build/opengl_simulation`: run the desktop app on Linux/macOS.
- `cmake -S . -B build -DCMAKE_TOOLCHAIN_FILE=/path/to/vcpkg/scripts/buildsystems/vcpkg.cmake`: reuse an existing `vcpkg` install instead of the local bootstrap.

## Coding Style & Naming Conventions
Use C++17. Follow the existing style in `src/`: 4-space indentation, opening braces on the same line for functions and control blocks, and descriptive comments only where the math or rendering flow is non-obvious. Keep type names in `PascalCase` (`Vector2`), functions in `camelCase` (`framebuffer_size_callback` is an existing exception from GLFW callback naming), and constants/global state readable and explicit. Prefer small, mode-specific functions over adding more logic to `main.cpp`.

## Testing Guidelines
There is no automated test suite yet. Treat a clean configure/build as the minimum gate, then run the app and sanity-check all three modes: Reflection, Refraction, and Prism. When adding non-trivial math, verify angle calculations and ray paths interactively. If you introduce tests later, place them in a top-level `tests/` directory and wire them into CMake.

## Commit & Pull Request Guidelines
History is minimal and uses imperative subjects, for example: `Implement ray tracing simulation with reflection and refraction`. Keep commit messages short, imperative, and feature-focused. For pull requests, include:
- a brief summary of behavior changes
- build/run confirmation on your platform
- screenshots or short recordings for UI/visual changes
- linked issues when applicable

## Dependency & Configuration Notes
Do not commit `.vcpkg/`, `build/`, or `vcpkg_installed/`. Update [`vcpkg.json`](./vcpkg.json) together with [`README.md`](./README.md) when dependency workflow changes.
