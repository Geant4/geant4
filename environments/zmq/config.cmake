#==============================================================================
# CMake configuration
#
# NOTE:
# CACHE variables can be changed in CMake CLI with -D option.
#==============================================================================

set(CMAKE_INSTALL_PREFIX /zzz
    CACHE STRING "Install prefix")

# Geant4 installation path
set(GEANT4_INSTALL /zzz
    CACHE STRING "Geant4 installation path")

# visualization flag
set(ENABLE_VIS FALSE CACHE BOOL "Enable visualization flag")

# Optimizaton / Debug flags
set(OPTIMIZE TRUE CACHE BOOL "Optimizaton flag (O3)")
set(DEBUG FALSE CACHE BOOL "Debug mode")

# Development flag (set false for release)
set(DEVMODE FALSE CACHE BOOL "Development mode")

# Optional configurations
# Compiler
#set(CMAKE_C_COMPILER /opt/gcc-4.8.5/bin/gcc)
#set(CMAKE_CXX_COMPILER /opt/gcc-4.8.5/bin/g++)
