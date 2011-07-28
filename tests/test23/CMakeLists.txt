cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test23)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test23 test23.cc src/*.cc)
GEANT4_EXECUTABLE(test23_hadronic_exerciser test23_hadronic_exerciser.cc)

#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test23_hadronic_exerciser COMMAND test23_hadronic_exerciser
                                          OUTPUT  test23-had.in
                                          BUILD test23_hadronic_exerciser)
GEANT4_ADD_TEST(test23 COMMAND test23 ${CMAKE_CURRENT_SOURCE_DIR}/test23.in
                       BUILD test23
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(test23-had COMMAND test23 test23-had.in
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test23_hadronic_exerciser)


