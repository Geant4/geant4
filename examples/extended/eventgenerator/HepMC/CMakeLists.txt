#---Adding all HepMC examples subdirectories explicitly 

cmake_minimum_required(VERSION 3.16...3.21)

find_package(Geant4)
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# HepMC examples require HepMC
#
find_package(HepMC REQUIRED)

add_subdirectory(HepMCEx01)
add_subdirectory(HepMCEx02)
add_subdirectory(MCTruth)
