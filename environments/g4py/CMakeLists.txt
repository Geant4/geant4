# - Top-level CMakeLists.txt for Geant4Py

cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

#------------------------------------------------------------------------------
project(Geant4Py)
#------------------------------------------------------------------------------
# installation prefixes for libraries
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR})

# debug mode
set(DEBUG FALSE CACHE BOOL "Debug Mode (Debug On)")

enable_testing()

#------------------------------------------------------------------------------
# Do not edit below
#------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 11)

# cmake modules
set(CMAKE_MODULE_PATH
    ${PROJECT_SOURCE_DIR}/cmake/Modules
    ${CMAKE_MODULE_PATH})

include(SetInstallPath)

# searching packages...
find_package(Geant4 REQUIRED)
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(Boost)
find_package(XercesC)
find_package(ROOT)

if(NOT GEANT4_FOUND)
  message(FATAL_ERROR "NOT Found Geant4.")
endif()

#------------------------------------------------------------------------------
# parameters for building
message(STATUS "--------------------------------------------------------")
message(STATUS "Parameters for building")

# build options
if(NOT CMAKE_BUILD_TYPE)
  if(DEBUG)
    set(CMAKE_BUILD_TYPE "Debug")
  else()
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
  endif()
endif()
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")

message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "--------------------------------------------------------")

#------------------------------------------------------------------------------
# add subdiretories...

# adding sub directories to the project
# libraries
add_subdirectory(source)
add_subdirectory(site-modules)

# examples
add_subdirectory(examples)

# tests
add_subdirectory(tests)
