#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(AnaEx01)

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
# Locate sources and headers for this project
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include 
                    # uncomment the line below to include HBOOK analysis
                    #${CMAKE_CURRENT_SOURCE_DIR}/../../common/analysis/include
                    ${CMAKE_CURRENT_SOURCE_DIR}/../shared/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

#----------------------------------------------------------------------------
# Get examples sources from common and shared
#
list(APPEND sources 
  # uncomment the line below to include HBOOK analysis
  # ${CMAKE_CURRENT_SOURCE_DIR}/../../common/analysis/src/ExG4HbookAnalysisManager.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/../shared/src/DetectorConstruction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/../shared/src/DetectorMessenger.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/../shared/src/EventAction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/../shared/src/PhysicsList.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/../shared/src/PrimaryGeneratorAction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/../shared/src/RunAction.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/../shared/src/SteppingAction.cc)
  
#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(AnaEx01 EXCLUDE_FROM_ALL AnaEx01.cc ${sources})
#target_link_libraries(AnaEx01 ${Geant4_LIBRARIES} -lXt)
target_link_libraries(AnaEx01 ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(ANAEX01_SCRIPTS
  AnaEx01.in
  #AnaEx01.out
  #init.mac
  #init_vis.mac
  run.mac
  vis.mac
  )

foreach(_script ${ANAEX01_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()


