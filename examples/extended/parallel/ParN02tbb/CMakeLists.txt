# - CMakeLists.txt for Geant4 MT Example ParN02
#
#-----------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ParN02tbb)
add_definitions("-g -O0 -DTBBDEBUG")
#-----------------------------------------------------------------------
# Find Geant4MT package, using the 'mt' component to enforce this.
# As with other examples, we activate all available UI and Vis drivers by 
# default.
# You can set WITH_GEANT4_UIVIS to OFF via the command line or 
# ccmake/cmake-gui to build a batch mode only executable
#
#option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)

find_package(Geant4 REQUIRED)

find_package(TBB REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${TBB_INCLUDE_DIRS} )
message("ParN02tbb : Using include " ${TBB_INCLUDE_DIRS} " for TBB" )

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so that they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(ParN02tbb ParN02tbb.cc ${sources} ${headers})
target_link_libraries(ParN02tbb ${Geant4_LIBRARIES} ${TBB_LIBRARIES} ${TBB_MALLOC_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build ParN02tbb. This is so that we can run the executable directly 
# because it relies on these scripts being in the current working 
# directory.
#
set(PARN02TBB_SCRIPTS
  run1debug.mac
  run1.mac
  )

foreach(_script ${PARN02TBB_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(N02tbb DEPENDS ParN02tbb)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ParN02tbb DESTINATION bin)

