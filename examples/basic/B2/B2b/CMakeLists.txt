cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(B2b)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

add_executable(exampleB2b exampleB2b.cc ${sources})
target_link_libraries(exampleB2b ${Geant4_LIBRARIES})

install(TARGETS exampleB2b DESTINATION bin )

