#---Adding all HepMC3 examples subdirectories explicitly 

cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})
find_package(Geant4)
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# HepMC examples require HepMC
#
find_package(HepMC3 REQUIRED)

add_subdirectory(HepMCEx01)
add_subdirectory(HepMCEx02)
add_subdirectory(MCTruth)
