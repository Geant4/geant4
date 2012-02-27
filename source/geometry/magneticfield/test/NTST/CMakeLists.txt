cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(testNTST)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

# Template
GEANT4_EXECUTABLE(testNTST testNTST src/*.cc)
#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(testNTST-2xa COMMAND testNTST ${CMAKE_CURRENT_SOURCE_DIR}/run2xa.mac
                       BUILD testNTST)
GEANT4_ADD_TEST(testNTST-2xb COMMAND testNTST ${CMAKE_CURRENT_SOURCE_DIR}/run2xb.mac
                       BUILD testNTST)
GEANT4_ADD_TEST(testNTST-2xc COMMAND testNTST ${CMAKE_CURRENT_SOURCE_DIR}/run2xb.mac
                       BUILD testNTST)
