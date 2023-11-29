#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(${common_prefix}common)

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
# Option to build examples in addition to test
#
option(COMMON_EXAMPLES "Activate building example executables" ON)

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
# Add the library
#
add_library(${common_prefix}common ${sources} ${headers})
target_link_libraries(${common_prefix}common ${Geant4_LIBRARIES} )

#----------------------------------------------------------------------------
# Add the test executable, and link it to the Geant4 libraries and common library
#
add_executable(${common_prefix}testCommon testCommon.cc )
target_link_libraries(${common_prefix}testCommon ${common_prefix}common ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Optionally, add the example executables, and link it to the Geant4 libraries and common library
#
if(COMMON_EXAMPLES)
  add_executable(exampleCommon exampleCommon.cc )
  target_link_libraries(exampleCommon ${common_prefix}common ${Geant4_LIBRARIES})
  add_executable(exampleCommonNoVis exampleCommonNoVis.cc )
  target_link_libraries(exampleCommonNoVis ${common_prefix}common ${Geant4_LIBRARIES})
endif(COMMON_EXAMPLES)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build detectorConstruction. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(detectorConstruction_SCRIPTS
  init.mac init_vis.mac run.mac vis.mac
)

foreach(_script ${detectorConstruction_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
  )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
# and the library to 'lib' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ${common_prefix}testCommon DESTINATION bin)
if (COMMON_EXAMPLES)
  install(TARGETS exampleCommon DESTINATION bin)
  install(TARGETS exampleCommonNoVis DESTINATION bin)
endif(COMMON_EXAMPLES)
install(TARGETS ${common_prefix}common DESTINATION lib)
