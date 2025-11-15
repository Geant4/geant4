# Build System Reference

## Overview

Geant4 uses CMake as its build system, providing cross-platform support and flexible configuration options.

## CMake Requirements

**Minimum Version**: 3.16
**Maximum Tested Version**: 3.27

From CMakeLists.txt:19:
```cmake
cmake_minimum_required(VERSION 3.16...3.27)
```

## Project Configuration

### Project Definition

The project is defined in CMakeLists.txt:29-35:

```cmake
project(Geant4
  DESCRIPTION "C++ toolkit for simulating the passage of particles through matter"
  HOMEPAGE_URL "https://geant4.cern.ch")
set(${PROJECT_NAME}_VERSION_MAJOR 11)
set(${PROJECT_NAME}_VERSION_MINOR  4)
set(${PROJECT_NAME}_VERSION_PATCH  0)
set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}")
```

### Version Information

- **Major**: 11
- **Minor**: 4
- **Patch**: 0
- **Full Version**: 11.4.0

## Build Requirements

### Out-of-Source Builds

Geant4 **requires** out-of-source builds. In-source builds are explicitly forbidden.

From CMakeLists.txt:6-14:

```cmake
if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  message(STATUS "Geant4 requires an out-of-source build.")
  message(STATUS "Please remove these files from ${CMAKE_CURRENT_BINARY_DIR} first:")
  message(STATUS "CMakeCache.txt")
  message(STATUS "CMakeFiles")
  message(STATUS "Once these files are removed, create a separate directory")
  message(STATUS "and run CMake from there")
  message(FATAL_ERROR "in-source build detected")
endif()
```

### Correct Build Process

```bash
# Create a build directory
mkdir build
cd build

# Configure
cmake ..

# Build
make -j$(nproc)

# Install (optional)
make install
```

### Incorrect Build Process

```bash
# This will FAIL - in-source build not allowed
cmake .
```

## Custom CMake Modules

Geant4 provides custom CMake modules in `cmake/Modules/`:

From CMakeLists.txt:41-43:
```cmake
set(CMAKE_MODULE_PATH
  ${PROJECT_SOURCE_DIR}/cmake/Modules
  ${CMAKE_MODULE_PATH})
```

### Custom Make Rules

Custom C++ compiler rules are defined in:
- `cmake/Modules/G4MakeRules_cxx.cmake`

From CMakeLists.txt:22:
```cmake
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/G4MakeRules_cxx.cmake)
```

## Main CMake Module

The main CMake configuration is included from:

CMakeLists.txt:49:
```cmake
include(G4CMakeMain)
```

This module contains the detailed build logic and is located in `cmake/Modules/G4CMakeMain.cmake`.

## Build Directories

### Source Build

The main source code build configuration is in:
- `source/CMakeLists.txt`

### Examples

Example applications can be built separately and are located in:
- `examples/`

## Configuration Options

CMake configuration options can be passed during the configuration step:

```bash
cmake .. -DOPTION_NAME=VALUE
```

For detailed configuration options, refer to the [Installation Guide](http://cern.ch/geant4-userdoc/UsersGuides/InstallationGuide/html).

## Platform Support

Geant4 supports:
- Linux
- macOS
- Windows

The CMake build system handles platform-specific configuration automatically.

## Related Files

| File | Purpose |
|------|---------|
| `CMakeLists.txt` | Top-level build configuration |
| `cmake/Modules/` | Custom CMake modules |
| `source/CMakeLists.txt` | Source code build configuration |
| `config/` | Additional configuration files |

## See Also

- [Installation Guide](http://cern.ch/geant4-userdoc/UsersGuides/InstallationGuide/html)
- [Getting Started](../getting-started)
- [Architecture Overview](../architecture)
