cmake_minimum_required(VERSION 2.8)

#----------------------------------------------------------------------------

project(chem4)

#----------------------------------------------------------------------------
# Build plot program if ROOT built with GUI via CMake is found
#
find_package(ROOT COMPONENTS Gui QUIET)
#message(STATUS ROOT_USE_FILE " ${ROOT_USE_FILE}")
if(ROOT_FOUND AND ROOT_USE_FILE)
  set(BUILD_CHEM4_PLOT ON)
  add_subdirectory(plot)
else()
  message("-- G4 Examples : The chem4/plot requires ROOT built with GUI via CMake, plot is disabled")
endif()

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
add_executable(chem4 chem4.cc ${sources} ${headers})
target_link_libraries(chem4 ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build chem4_proj. This is so that we can run the executable directly because 
# it relies on these scripts being in the current working directory.
#
file(GLOB CHEM4_SCRIPTS
  ${PROJECT_SOURCE_DIR}/*.in
  ${PROJECT_SOURCE_DIR}/*.mac)
#message(STATUS CHEM4_SCRIPTS " ${CHEM4_SCRIPTS}")

foreach(_script ${CHEM4_SCRIPTS})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/.
    COPYONLY
    )
endforeach()

source_group("macros" FILES ${CHEM4_SCRIPTS})
target_sources(chem4 PUBLIC ${CHEM4_SCRIPTS})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS chem4 DESTINATION bin )

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
#project(chem4_proj)
#add_custom_target(chem4_proj DEPENDS chem4 plotG)

if(BUILD_CHEM4_PLOT)
  add_dependencies(plotG chem4)
endif()
