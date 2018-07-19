#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(neuron)

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
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(neuron neuron.cc ${sources} ${headers})
target_link_libraries(neuron ${Geant4_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build neuron. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
file(GLOB MAC_FILES ${PROJECT_SOURCE_DIR}/*.mac 
                    ${PROJECT_SOURCE_DIR}/*.in
                    ${PROJECT_SOURCE_DIR}/*.C
                    ${PROJECT_SOURCE_DIR}/README
                    ${PROJECT_SOURCE_DIR}/*.swc
                    ${PROJECT_SOURCE_DIR}/*.dat
                    )

foreach(_script ${MAC_FILES})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/.
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS neuron DESTINATION bin)

