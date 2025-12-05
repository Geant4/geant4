# ----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(microtrack)

# ----------------------------------------------------------------------------
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

# ----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

# ----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIRS})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# ----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(microtrack microtrack.cc ${sources} ${headers})
target_link_libraries(microtrack ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

# ----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we build
# microtrack. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(microtrack_SCRIPTS run.mac init_vis.mac vis.mac processMicrotrack.C)

foreach(_script ${microtrack_SCRIPTS})
  configure_file(${PROJECT_SOURCE_DIR}/${_script}
                 ${PROJECT_BINARY_DIR}/${_script} COPYONLY)
endforeach()

# ----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS microtrack DESTINATION bin)
