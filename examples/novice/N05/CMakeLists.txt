cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(N05)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

add_executable(exampleN05 EXCLUDE_FROM_ALL exampleN05.cc ${sources})
target_link_libraries(exampleN05 ${Geant4_LIBRARIES})

