#----------------------------------------------------------------------------
# Setup the project

cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(saxs)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable

option(WITH_GEANT4_UIVIS "Build SAXS with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries

add_executable(saxs saxs.cc ${sources} ${headers})
target_link_libraries(saxs ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build SAXS. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.

set(SAXS_SCRIPTS
  	init_vis.mac
  	vis.mac
	gui.mac
	geom.mac
	myFF.dat
	saxs.in
	saxs_slits.in
	test.in
	scattAnalysis.C
	ADXRD.C
)

foreach(_script ${SAXS_SCRIPTS})
	configure_file(
    	${PROJECT_SOURCE_DIR}/${_script}
    	${PROJECT_BINARY_DIR}/${_script}
    	COPYONLY
	)
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
# add_custom_target(saxs_custom DEPENDS saxs)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX

install(TARGETS saxs DESTINATION bin)

