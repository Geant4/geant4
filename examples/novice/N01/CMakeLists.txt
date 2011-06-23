cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(N01)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_LINK_LIBRARY(exampleN01 src/*.cc)
GEANT4_EXECUTABLE(exampleN01 exampleN01.cc LIBRARIES exampleN01)





