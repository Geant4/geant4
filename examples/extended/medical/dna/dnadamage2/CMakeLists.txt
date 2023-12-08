cmake_minimum_required(VERSION 3.16...3.27)

#----------------------------------------------------------------------------

project(dnadamage2)

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
include_directories(${PROJECT_SOURCE_DIR}/include)

include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)

file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(dnadamage2 dnadamage2.cc ${sources} ${headers})
target_link_libraries(dnadamage2 ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build dnadamage2_proj. This is so that we can run the executable directly because 
# it relies on these scripts being in the current working directory.
#
file(GLOB dnadamage2_SCRIPTS
  ${PROJECT_SOURCE_DIR}/*.in
  ${PROJECT_SOURCE_DIR}/*.mac
  ${PROJECT_SOURCE_DIR}/*.xyz
  ${PROJECT_SOURCE_DIR}/*.txt
  ${PROJECT_SOURCE_DIR}/*.src)
#message(STATUS dnadamage2_SCRIPTS " ${dnadamage2_SCRIPTS}")

foreach(_script ${dnadamage2_SCRIPTS})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/.
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS dnadamage2 DESTINATION bin )

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
project(dnadamage2_proj)
add_custom_target(dnadamage2_proj DEPENDS dnadamage2)
