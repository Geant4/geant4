#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(AnaEx03)

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
# Find AIDA (required package)
#
find_package(AIDA QUIET)
if(NOT AIDA_FOUND)
  message(STATUS "G4 Examples: AIDA package not found. --> AnaEx03 example disabled")  
  return()
else()
  add_definitions(-DG4ANALYSIS_USE) 
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/../shared/include 
                    ${Geant4_INCLUDE_DIR}
                    ${AIDA_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Get examples sources from common and shared
#
list(APPEND sources 
  ${PROJECT_SOURCE_DIR}/../shared/src/DetectorConstruction.cc
  ${PROJECT_SOURCE_DIR}/../shared/src/DetectorMessenger.cc
  ${PROJECT_SOURCE_DIR}/../shared/src/EventAction.cc
  ${PROJECT_SOURCE_DIR}/../shared/src/PhysicsList.cc
  ${PROJECT_SOURCE_DIR}/../shared/src/PrimaryGeneratorAction.cc
  ${PROJECT_SOURCE_DIR}/../shared/src/RunAction.cc
  ${PROJECT_SOURCE_DIR}/../shared/src/SteppingAction.cc)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(AnaEx03 AnaEx03.cc ${sources} ${headers})
target_link_libraries(AnaEx03 ${Geant4_LIBRARIES} ${AIDA_LIBRARIES} )

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build AnaEx03. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(AnaEx03_SCRIPTS
    AnaEx03.in run.mac vis.mac
  )

foreach(_script ${AnaEx03_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS AnaEx03 DESTINATION bin)

