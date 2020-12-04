#---Adding all HepMC examples subdirectories explicitly 

cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

find_package(Geant4)
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# HepMC examples require HepMC
#
find_package(HepMC REQUIRED)

add_subdirectory(HepMCEx01)
add_subdirectory(HepMCEx02)
add_subdirectory(MCTruth)
