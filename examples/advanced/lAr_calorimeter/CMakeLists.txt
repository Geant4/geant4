#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(lArCal)

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

find_package(AIDA QUIET)
if(NOT AIDA_FOUND)
  message(STATUS "G4 Examples: AIDA package not found. --> lArCalorimeter example disabled")  
  return()
else()
  add_definitions(-DG4ANALYSIS_USE) 
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Setup include directory for this project
include_directories(${PROJECT_SOURCE_DIR}/include ${AIDA_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(lArCal lArCal.cc ${sources} ${headers})
target_link_libraries(lArCal ${Geant4_LIBRARIES} ${AIDA_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build lArCal. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLELARCALORIMETER_SCRIPTS
  gui.mac
  run1.mac
  run2.mac
  vis.mac
  icons.mac
  )

foreach(_script ${EXAMPLELARCALORIMETER_SCRIPTS})
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
add_custom_target(ExamplelArCal DEPENDS lArCal)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS lArCal DESTINATION bin)


