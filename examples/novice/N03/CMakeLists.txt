#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(N03)

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
add_executable(exampleN03 exampleN03.cc ${sources} ${headers})
target_link_libraries(exampleN03 ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build N03. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEN03_SCRIPTS
  exampleN03.err
  exampleN03.in
  exampleN03.large_N.in
  exampleN03.out
  exo1.mac
  exo2.mac
  gui.mac
  newgeom.mac
  rndmSeed.mac
  run1.mac
  run2.mac
  vis.mac
  visTutor/exN03Tree0.mac   
  visTutor/exN03Vis13.loop  
  visTutor/heprep.mac
  visTutor/exN03Vis0.mac    
  visTutor/exN03Vis13.mac   
  visTutor/History
  visTutor/exN03Vis12.loop  
  visTutor/exN03Vis5.mac
  visTutor/exN03Vis12.mac   
  visTutor/g4Xt.xrm
  )

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/visTutor)

foreach(_script ${EXAMPLEN03_SCRIPTS})
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
add_custom_target(N03 DEPENDS exampleN03)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS exampleN03 DESTINATION bin)


