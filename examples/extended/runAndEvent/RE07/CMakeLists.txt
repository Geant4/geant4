#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(RE07)

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
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Define sources for this project
#
set(sources
  exampleRE07.cc
  src/ActionInitialization.cc
  src/DetectorConstruction.cc
  src/DetectorMessenger.cc
  src/EmStandardPhysicsTrackingManager.cc
  src/EventAction.cc
  src/PhysicsList.cc
  src/PhysicsListEmSpecialized.cc
  src/PhysicsListEmStandardTracking.cc
  src/PhysicsListMessenger.cc
  src/PrimaryGeneratorAction.cc
  src/RunAction.cc
  src/Run.cc
  src/SpecializedTrackingManager.cc
  src/SteppingAction.cc
  src/TrackingAction.cc)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exampleRE07 ${sources})
target_link_libraries(exampleRE07 ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS exampleRE07 DESTINATION bin)
