#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(STCyclotron)
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
include_directories(${PROJECT_SOURCE_DIR}/include)
#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})
#---------------------------------------------------
# Locate sources and headers for this project
#     
include_directories(${PROJECT_SOURCE_DIR}/include 
  ${Geant4_INCLUDE_DIR})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(STCyclotron STCyclotron.cc ${sources} ${headers})

target_link_libraries(STCyclotron ${Geant4_LIBRARIES}  ${ROOT_LIBRARIES})
#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build rdecay02. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(STCYCLOTRON_SCRIPTS
   init.mac
   init_parameters.mac
   run.mac
   vis.mac
   Macro/Vis/vis.mac
   Macro/Vis/init_vis.mac
   Macro/init.mac
   Macro/init_parameters.mac
   Macro/run.mac
   Macro/Vis/vis.mac
   Macro/Vis/init_vis.mac
   Plot.C
   Macro/GUI/gui.mac
   Macro/GUI/icons.mac
   Macro/GUI/run.png
   Macro/Material/Target/Ni64_enriched_95.mac
   Macro/Material/Target/Ni64_pure.mac
   Macro/Material/Target/My_designed_material.mac
   )

foreach(_script ${STCYCLOTRON_SCRIPTS})
 configure_file(
   ${PROJECT_SOURCE_DIR}/${_script}
   ${PROJECT_BINARY_DIR}/${_script}
   COPYONLY
   )
endforeach()

