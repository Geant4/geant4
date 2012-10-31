 #----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(name DMX)
project(${name})

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#Check here for DMXENV_GPS_USE
#Use GPS if the environment variable is set, and the user has *not* supplied the CMake 
#command line argument. Otherwise, the argument of the command line has the priority.
#
if("$ENV{DMXENV_GPS_USE}" AND NOT DEFINED DMXENV_GPS_USE)
  set(DMXENV_GPS_USE_DEFAULT ON)
else()
  set(DMXENV_GPS_USE_DEFAULT ${DMXENV_GPS_USE})
endif()

option(DMXENV_GPS_USE "Build DMX with support for GPS" ${DMXENV_GPS_USE_DEFAULT})
mark_as_advanced(DMXENV_GPS_USE)

if(DMXENV_GPS_USE)
  add_definitions(-DDMXENV_GPS_USE)
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
# Setup include directory for this project
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(${name} EXCLUDE_FROM_ALL ${name}.cc ${sources})
target_link_libraries(${name} ${Geant4_LIBRARIES})


