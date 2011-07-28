cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test202)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test202 test202.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test202-build BUILD test202)
