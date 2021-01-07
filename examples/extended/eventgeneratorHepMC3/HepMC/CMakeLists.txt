#---Adding all HepMC examples subdirectories explicitly 

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})
find_package(Geant4)
include(${Geant4_USE_FILE})
set ( USE_HEPMC3 FALSE)
set ( USE_HEPMC3 TRUE)
add_subdirectory(HepMCEx01)
add_subdirectory(HepMCEx02)
add_subdirectory(MCTruth)
