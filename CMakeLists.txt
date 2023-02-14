#-----------------------------------------------------------------------
# - Top Level CMakeLists.txt for Geant4 Build
#-----------------------------------------------------------------------
# - Enforce an out-of-source builds before anything else
#
if(${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  message(STATUS "Geant4 requires an out-of-source build.")
  message(STATUS "Please remove these files from ${CMAKE_CURRENT_BINARY_DIR} first:")
  message(STATUS "CMakeCache.txt")
  message(STATUS "CMakeFiles")
  message(STATUS "Once these files are removed, create a separate directory")
  message(STATUS "and run CMake from there")
  message(FATAL_ERROR "in-source build detected")
endif()

#-----------------------------------------------------------------------
# - Define CMake requirements and override make rules as needed
#
cmake_minimum_required(VERSION 3.16...3.24)

# - Make overrides for default flags, so they appear in interfaces
set(CMAKE_USER_MAKE_RULES_OVERRIDE_CXX ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/G4MakeRules_cxx.cmake)

#-----------------------------------------------------------------------
# - Project definition and basic configuration
#   Version handled manually as project(... VERSION ...) is not used
#   in tests/examples which are subprojects. All calls must use
#   the same form.
project(Geant4
  DESCRIPTION "C++ toolkit for simulating the passage of particles through matter"
  HOMEPAGE_URL "https://geant4.cern.ch")
set(${PROJECT_NAME}_VERSION_MAJOR 11)
set(${PROJECT_NAME}_VERSION_MINOR  1)
set(${PROJECT_NAME}_VERSION_PATCH  1)
set(${PROJECT_NAME}_VERSION "${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH}")

# - Prepend our own CMake Modules to the search path
# NB: if our custom modules include others that we don't supply, those in
# the base path will be used, so watch for incompatibilities!!
#
set(CMAKE_MODULE_PATH
  ${PROJECT_SOURCE_DIR}/cmake/Modules
  ${CMAKE_MODULE_PATH})

#-----------------------------------------------------------------------
# - Include CMake category main module
#   Factored into category for convenience in tagging
#
include(G4CMakeMain)

