#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.16...3.21)
project(ICRP145phantoms)

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
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(ICRP145phantoms ICRP145phantoms.cc ${sources} ${headers})
target_link_libraries(ICRP145phantoms ${Geant4_LIBRARIES})

# Depend on data for runtime
add_dependencies(ICRP145phantoms ICRP145data)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory. This is so that we can run the 
# executable directly because it relies on these scripts being in the current 
# working directory.
#
set(EXTERNAL_SCRIPTS
  example.in
  init_vis.mac
  source.mac
  vis.mac
  )

foreach(_script ${EXTERNAL_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

include(ExternalProject)
ExternalProject_Add(ICRP145data
  SOURCE_DIR ${PROJECT_BINARY_DIR}/ICRP145data
  URL https://cern.ch/geant4-data/datasets/examples/advanced/ICRP145Phantoms/ICRP145data.tar.gz
  URL_MD5 5d2a46eaa416349ff7d0f3415345f886
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  )

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ICRP145phantoms DESTINATION bin)
