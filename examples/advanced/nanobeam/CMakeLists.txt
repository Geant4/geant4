#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(nanobeam)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI drivers by default
# You can set WITH_GEANT4_UI to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UI "Build example with Geant4 UI drivers" ON)
if(WITH_GEANT4_UI)
  find_package(Geant4 REQUIRED ui_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

if(Geant4_builtin_clhep_FOUND)
  message(STATUS "G4 Examples: CLHEP external package not found  --> nanobeam disabled")
  return()
endif() 
#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(nanobeam nanobeam.cc ${sources} ${headers})
target_link_libraries(nanobeam ${Geant4_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build nanobeam. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(nanobeam_SCRIPTS
    default.mac nanobeam.out plot.C
  )

foreach(_script ${nanobeam_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS nanobeam DESTINATION bin)

