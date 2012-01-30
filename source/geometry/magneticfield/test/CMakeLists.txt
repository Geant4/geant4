cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(test-magfield)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_subdirectory(NTST)
add_subdirectory(field07)

# Template
# GEANT4_EXECUTABLE(test18 test18.cc src/*.cc)

GEANT4_EXECUTABLE(testPropagateMagField  testPropagateMagField.cc)
GEANT4_EXECUTABLE(testProElectroMagField testProElectroMagField.cc)
GEANT4_EXECUTABLE(testProPerpSpin testProPerpSpin.cc)
GEANT4_EXECUTABLE(testPropagateSpin testPropagateSpin.cc)
