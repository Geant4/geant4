cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(B3)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

add_executable(exampleB3 exampleB3.cc ${sources})
target_link_libraries(exampleB3 ${Geant4_LIBRARIES})

add_custom_target(B3 DEPENDS exampleB3)

install(TARGETS exampleB3 DESTINATION bin )
