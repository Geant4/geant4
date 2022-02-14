#----------------------------------------------------------------------------
# Setup the project

cmake_minimum_required(VERSION 3.16...3.21)

project(py8decayer)

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
# Find Pythia8 (required package)
#
# Ideally, should be done via find_package, e.g.
find_package(Pythia8 REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
		    ${PYTHIA8_INCLUDES})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/src/*.c)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(pythia8_decayer pythia8_decayer.cc ${sources} ${headers})
target_link_libraries(pythia8_decayer  ${Geant4_LIBRARIES} ${PYTHIA8_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build decayer6. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
# ---> NO NEED SO FAR !!!
#
#set(py8decayer_SCRIPTS
#    init_vis.mac pythia8_decayer.in pythia8_decayer.out vis.mac
#  )

#foreach(_script ${py8decayer_SCRIPTS})
#  configure_file(
#    ${PROJECT_SOURCE_DIR}/${_script}
#    ${PROJECT_BINARY_DIR}/${_script}
#    COPYONLY
#    )
#endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS pythia8_decayer DESTINATION bin)

