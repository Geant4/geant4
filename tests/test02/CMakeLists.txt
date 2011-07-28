cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test02)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test02 test02.cc src/*.cc)
GEANT4_EXECUTABLE(test02_had test02_hadronic_exerciser.cc)


#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test02 COMMAND test02 ${CMAKE_CURRENT_SOURCE_DIR}/test02.in
                       BUILD test02
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(test02-largeN COMMAND test02 ${CMAKE_CURRENT_SOURCE_DIR}/test02.large_N.in
                       DEPENDS test02
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

GEANT4_ADD_TEST(test02-gen COMMAND test02_had OUTPUT test02-had.in
                           BUILD test02_had)
GEANT4_ADD_TEST(test02-had COMMAND test02 test02-had.in
                           DEPENDS test02-gen test02
                           ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
