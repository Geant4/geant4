#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(Hadr02)

#----------------------------------------------------------------------------
# Prepend our own CMake Modules to the search path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})

#----------------------------------------------------------------------------
# Check whether CRMC should be used or not

option(G4_USE_CRMC "Using CRMC" OFF) 
if(G4_USE_CRMC)
  message(STATUS "G4_USE_CRMC ON : Using CRMC for building Hadr02")
  add_definitions(-DG4_USE_CRMC)
  message(STATUS "CMAKE_MODULE_PATH =" ${CMAKE_MODULE_PATH} " ; PROJECT_SOURCE_DIR = " ${PROJECT_SOURCE_DIR})
  find_package(CRMC REQUIRED)
  if(CRMC_FOUND)
    message(STATUS "CRMC was found :" ${CRMC_INCLUDE_DIR})
  else()
    message(FATAL_ERROR "CRMC was NOT found !")
  endif()
else()
  message(STATUS "G4_USE_CRMC OFF : NOT using CRMC for building Hadr02")
endif()

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

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${CRMC_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${CRMC_INCLUDE_DIR}/*.h ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(Hadr02 Hadr02.cc ${sources} ${headers})
target_link_libraries(Hadr02 ${CRMC_LIBRARIES} ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build Hadr02. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(Hadr02_SCRIPTS
    hadr02.in hadr02.out urqmd.in crmc.in vis.mac
  )

foreach(_script ${Hadr02_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

if(G4_USE_CRMC)
  if(CRMC_FOUND)
    configure_file(${CRMC_BUILD_DIR}/crmc.param ${PROJECT_BINARY_DIR} COPYONLY)
  endif()
endif()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS Hadr02 DESTINATION bin)

