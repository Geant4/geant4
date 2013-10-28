#----------------------------------------------------------------------------
#
#  Note-1: Env.variables G4INSTALL, G4LIB, ROOTSYS need to be setup
#
#  Note-2: Maybe ROOT can be pulled in via ROOT_INCLUDE_DIR - need to double check...
#
# Usage:
#
# mkdir /tests/t23-bld
# cd /tests/t23-bld
# cmake -DGeant4_DIR=$G4LIB/Geant4-9.6.0 \
#       ../test23
# make
#
# Note-3: in some cases, one would also use -DGeant4_SOURCE_DIR=$G4INSTALL option,
#         but here's it's not necessary - the system will even tell about that
#


# Specify mandatory
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

#----Standalone--------------------------------------------------------------------------

#------------------------------------------------------------------------------
# - Enforce an out-of-source builds before anything else
#
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "Geant4 requires an out-of-source build.")
    message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
    message(STATUS "CMakeCache.txt")
    message(STATUS "CMakeFiles")
    message(STATUS "Once these files are removed, create a separate directory")
    message(STATUS "and run CMake from there")
    message(FATAL_ERROR "in-source build detected")
endif()

project(test23)
# message("Project Name: ${CMAKE_PROJECT_NAME}")

#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#
# find additional ROOT package
#
find_package(ROOT QUIET)
if(NOT ROOT_FOUND)
  message(STATUS "G4 TESTS: ROOT package not found. ")  
  return()
endif()

#
# an attempt to factor out and re-use common cmake fragment
# ... so far, doesn't work:
#     the inclusion is fine, and inside the fragment, things seem to startout fine
#     but I need to figure out how to properly return some ENV variables to the master
#
#include(${PROJECT_SOURCE_DIR}/../test23/shared/HeadersAndSources.cmake)
#HEADERS_AND_SOURCES( sources, headers )
#message("sources: ${sources}")

#----------------------------------------------------------------------------
# setup include directories
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/../test23/shared/g4app/include
                    ${Geant4_INCLUDE_DIR})

#----------------------------------------------------------------------------
#
# sources and headers for this project
#
file(GLOB prj_sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB prj_headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#
#  sources from shared code
#
file(GLOB g4app_sources 
  ${PROJECT_SOURCE_DIR}/../test23/shared/g4app/src/*.cc)
file(GLOB g4app_headers 
  ${PROJECT_SOURCE_DIR}/../test23/shared/g4app/include/*.hh)

#
# additional headers and spources from ROOT (if any)
#
include_directories(${PROJECT_SOURCE_DIR}/../test23/shared/rootanalysis/include
		    ${ROOTSYS}/include)
file(GLOB rootana_sources 
  ${PROJECT_SOURCE_DIR}/../test23/shared/rootanalysis/src/*.cc)
file(GLOB rootana_headers 
  ${PROJECT_SOURCE_DIR}/../test23/shared/rootana/include/*.hh)

#
# finally, put all headers and sources together 
#
list(APPEND sources ${prj_sources} ${g4app_sources} ${rootana_sources})
list(APPEND headers ${prj_headers} ${g4app_headers} ${rootana_headers})

#----------------------------------------------------------------------------
#
# specify what executable, and what are the sources for it...
#
add_executable(test23 test23.cc ${sources} ${headers})

#
# ... plus, the external libraries
#
target_link_libraries(test23 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

# copy over config for the ctest
#
set(test23_SCRIPTS ctest23-harp.config)
foreach(_script ${test23_SCRIPTS})
   configure_file(${PROJECT_SOURCE_DIR}/${_script} ${PROJECT_BINARY_DIR}/${_script} COPYONLY)
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS test23 DESTINATION bin)

