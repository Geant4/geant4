cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(exampleN03Con)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

add_executable(exampleN03Con EXCLUDE_FROM_ALL exampleN03Con.cc ${sources})
target_link_libraries(exampleN03Con ${Geant4_LIBRARIES})
