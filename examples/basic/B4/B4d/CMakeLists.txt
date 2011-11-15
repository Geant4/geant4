cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(B4d)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

add_executable(exampleB4d exampleB4d.cc ${sources})
target_link_libraries(exampleB4d ${Geant4_LIBRARIES})

install(TARGETS exampleB4d DESTINATION bin )

