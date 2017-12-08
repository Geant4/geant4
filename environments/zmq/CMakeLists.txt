#
#  CMakeLists.txt for buidlig zmq-geant4
#
#  Copyright 2017 Koichi Murakami
#
#  ***************************
#  !! DONOT EDIT THIS FILE  !!
#  ***************************
#
cmake_minimum_required(VERSION 3.3)
include(config.cmake)
#==============================================================================
project(ZMQ-GEANT4)
#==============================================================================

# version.h generation (true for release)
if(DEVMODE)
  set(SKIP_VERSION FALSE)
else()
  set(SKIP_VERSION TRUE)
endif()

# static / shared library
#set (LIBRARY_TYPE "STATIC")

# cmake modules
set(CMAKE_MODULE_PATH
    ${PROJECT_SOURCE_DIR}/cmake/Modules
    ${CMAKE_MODULE_PATH})

# library path (lib or lib64)
if(NOT DEFINED CMAKE_INSTALL_LIBDIR)
  set(_LIBDIR_DEFAULT "lib")
  if(CMAKE_SYSTEM_NAME MATCHES "Linux"
      AND NOT EXISTS "/etc/debian_version")
    if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
      set(_LIBDIR_DEFAULT "lib64")
    endif()
  endif()
  set(CMAKE_INSTALL_LIBDIR "${_LIBDIR_DEFAULT}")
endif()

# find packages
find_package(Geant4 REQUIRED)
find_package(ZMQ REQUIRED)

# c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# build options
if(NOT CMAKE_BUILD_TYPE)
  if(DEBUG)
    set(CMAKE_BUILD_TYPE "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG")
  elseif(OPTIMIZE)
    set(CMAKE_BUILD_TYPE "Release")
  else()
    set(CMAKE_BUILD_TYPE "RelWithDebInfo")
  endif()
endif()

#------------------------------------------------------------------------------
# parameters for building
message(STATUS "--------------------------------------------------------")
message(STATUS "Parameters for building")
message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "GEANT4_LIBRARY_DIR: ${GEANT4_LIBRARY_DIR}")

message(STATUS "--------------------------------------------------------")

#------------------------------------------------------------------------------
# version.h (#.##d# + rev# / dev flag d/a/b/rc)
#if(NOT SKIP_VERSION)
#set(G4ZMQ_VERSION_MAJOR 1)
#set(G4ZMQ_VERSION_MINOR 0.1)
#execute_process (COMMAND git rev-parse --short HEAD OUTPUT_VARIABLE
#                 G4ZMQ_VERSION_BUILD OUTPUT_STRIP_TRAILING_WHITESPACE)
#configure_file(
#  "${PROJECT_SOURCE_DIR}/cmake/version.h.in"
#  "${PROJECT_SOURCE_DIR}/zmq/version.h"
#)
#endif()

#------------------------------------------------------------------------------
# adding sub directories to the project
add_subdirectory(zmq)

# unit tests
if(DEVMODE)
include(CTest)
enable_testing()
#add_subdirectory(tests/ut)
endif()

# integration tests
#add_subdirectory(tests/it EXCLUDE_FROM_ALL)
