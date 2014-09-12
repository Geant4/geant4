# $Id$

# Common CMake configuration file for testAnalysis exacutables built
# with different compilation flags (TEST_ANALYSIS_XYZ).

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
#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B4c. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(TEST_ANALYSIS_SCRIPTS
  testAnalysis.in
  icons.mac
  gui.mac
  run.png
  init.mac
  init_vis.mac
  h1.mac
  h2.mac
  h3.mac
  p1.mac
  p2.mac
  run1.mac
  run2.mac
  vis.mac
  )

foreach(_script ${TEST_ANALYSIS_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/../macros/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

