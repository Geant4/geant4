#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(human_phantom)

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

#----------------------------------------------------------------------------
# Find gdml (optional package)
#
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
find_package(gdml QUIET)
if(GDML_FOUND)
  message(STATUS "G4 Examples: gdml found. --> human_phantom example with gdml enabled.") 
  # Uncomment this line if suitable (some customization may be needed)
  add_definitions(-DG4LIB_USE_GDML)   
else()
  message(STATUS "G4 Examples: gdml not found. --> human_phantom example with gdml disabled.")
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR}
                    ${GDML_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(phantom phantom.cc ${sources} ${headers})
target_link_libraries(phantom ${Geant4_LIBRARIES} ${GDML_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build human_phantom. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(human_phantom_SCRIPTS
    adultHead.mac adultMIRDFemale.mac adultMIRDMale.mac adultMIXFemale.mac adultORNLFemale.mac adultORNLMale.mac batch.mac dawnVis.mac default.mac noVis.mac openGLVis.mac run.mac vrmlVis.mac
  )

foreach(_script ${human_phantom_SCRIPTS})
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
add_custom_target(human_phantom DEPENDS phantom)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS phantom DESTINATION bin)

