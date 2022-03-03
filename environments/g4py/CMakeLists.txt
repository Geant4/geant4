# - Top-level CMakeLists.txt for Geant4Py
# Note that we are assuming it is always part of Geant4, never separately
# released. That's for later, if required at all.
# In particular, Geant4Py must be matched exactly to a given Geant4 release

# Require Python 3 or newer, calling Interp before Libs to help find a consistent set
# Boost python requirement is reliant on correct user config (not yet known how to check
# which Python an install of Boost Python uses)...
find_package(PythonInterp 3.0 REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(Boost 1.65 REQUIRED python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR})

# Variables and Functions to help build, organise, and install modules
include(${CMAKE_CURRENT_SOURCE_DIR}/G4PythonHelpers.cmake)

# - Build package/modules, tests if required
add_subdirectory(source)
if(GEANT4_ENABLE_TESTING)
  add_subdirectory(tests)
endif()

# Add the feature recording here so we can note the Python version
geant4_add_feature(GEANT4_USE_PYTHON "Building bindings for Python ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")

