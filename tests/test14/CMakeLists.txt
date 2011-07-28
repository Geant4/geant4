cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test14)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test14 test14.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test14 COMMAND test14 ${CMAKE_CURRENT_SOURCE_DIR}/test14.in
                       BUILD test14
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(test14-largeN COMMAND test14 ${CMAKE_CURRENT_SOURCE_DIR}/test14.large_N.in
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test14)
GEANT4_ADD_TEST(test14-lowE COMMAND test14 ${CMAKE_CURRENT_SOURCE_DIR}/test14_lowE.in
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test14)
GEANT4_ADD_TEST(test14-lowE-largeN COMMAND test14 ${CMAKE_CURRENT_SOURCE_DIR}/test14_lowE.large_N.in
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test14)
GEANT4_ADD_TEST(test14-penelope COMMAND test14 ${CMAKE_CURRENT_SOURCE_DIR}/test14_penelope.in
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test14)
GEANT4_ADD_TEST(test14-penelope-largeN COMMAND test14 ${CMAKE_CURRENT_SOURCE_DIR}/test14_penelope.large_N.in
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test14)
GEANT4_ADD_TEST(test14-polarised COMMAND test14 ${CMAKE_CURRENT_SOURCE_DIR}/test14_polarised.in
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test14)
GEANT4_ADD_TEST(test14-polarised-largeN COMMAND test14 ${CMAKE_CURRENT_SOURCE_DIR}/test14_polarised.large_N.in
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test14)
GEANT4_ADD_TEST(test14-std COMMAND test14 ${CMAKE_CURRENT_SOURCE_DIR}/test14_std.in
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test14)
GEANT4_ADD_TEST(test14-std-largeN COMMAND test14 ${CMAKE_CURRENT_SOURCE_DIR}/test14_std.large_N.in
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test14)


