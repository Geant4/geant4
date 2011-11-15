cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(B1)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

add_executable(exampleB1 exampleB1.cc ${sources})
target_link_libraries(exampleB1 ${Geant4_LIBRARIES})

add_custom_target(B1 DEPENDS exampleB1)

install(TARGETS exampleB1 DESTINATION bin )


