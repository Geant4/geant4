#---Adding all HepMC3 examples subdirectories explicitly 

cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})
find_package(Geant4)
include(${Geant4_USE_FILE})
add_subdirectory(HepMCEx01)
add_subdirectory(HepMCEx02)
add_subdirectory(MCTruth)
