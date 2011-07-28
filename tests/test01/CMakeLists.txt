cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test01)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test01 test01.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test01 COMMAND test01 ${CMAKE_CURRENT_SOURCE_DIR}/test01.in
                       BUILD test01)
GEANT4_ADD_TEST(test01-largeN COMMAND test01 ${CMAKE_CURRENT_SOURCE_DIR}/test01.large_N.in
                       DEPENDS test01)
