#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(CompositeCalorimeter)

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
add_executable(CompositeCalorimeter CompositeCalorimeter.cc ${sources} ${headers})
target_link_libraries(CompositeCalorimeter ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build CompositeCalorimeter. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLECOMPOSITECALORIMETER_SCRIPTS
  gui.mac
  test.g4mac
  vis.mac
  dataconf/g4testbeamhcal96.conf
  dataconf/testbeamhcal96.conf
  datageom/tbhcal96.geom
  datageom/tbhcal96hcal.geom
  datageom/tbhcal96xtal.geom
  dataglobal/fmap.tb96
  dataglobal/material.cms
  dataglobal/rotation.cms
  datavis/tbhcal96hcal.vis
  datavis/tbhcal96.vis
  datavis/tbhcal96xtal.vis
  )

foreach(_script ${EXAMPLECOMPOSITECALORIMETER_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(ExampleCompositeCalorimeter DEPENDS CompositeCalorimeter)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS CompositeCalorimeter DESTINATION bin)


