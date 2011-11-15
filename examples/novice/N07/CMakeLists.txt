cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(N07)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

add_executable(exampleN07 exampleN07.cc ${sources})
target_link_libraries(exampleN07 ${Geant4_LIBRARIES})

add_custom_target(N07 DEPENDS exampleN07)

install(TARGETS exampleN07 DESTINATION bin )

