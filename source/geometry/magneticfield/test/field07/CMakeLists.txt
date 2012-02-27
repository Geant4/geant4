cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(testNTST)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

# Template copied from test1x
## GEANT4_EXECUTABLE(field07 field07 src/*.cc)

# Template copied from exampleN02
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIR})
## file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh)
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)

add_executable(field07 EXCLUDE_FROM_ALL field07.cc ${sources})
target_link_libraries(field07 ${Geant4_LIBRARIES})

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(field07-run COMMAND field07 ${CMAKE_CURRENT_SOURCE_DIR}/field07.in
                       BUILD field07)
