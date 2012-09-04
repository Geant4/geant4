cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(mctruthex)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

find_package(HepMC QUIET)
if(NOT HEPMC_FOUND)
  message(STATUS "G4 Examples: HepMC package not found. --> MCTruth example disabled")  
  return()
endif()

include_directories(${PROJECT_SOURCE_DIR}/include ${HEPMC_INCLUDE_DIR})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(mctruthex mctruthex.cc ${sources} ${headers})
target_link_libraries(mctruthex ${Geant4_LIBRARIES} ${HEPMC_LIBRARIES})

