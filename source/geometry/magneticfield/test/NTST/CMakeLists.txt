cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(testNTST)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

# Template
GEANT4_EXECUTABLE(testNTST testNTST src/*.cc)
