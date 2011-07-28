cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test201)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test201 test201.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test201-build BUILD test201)
