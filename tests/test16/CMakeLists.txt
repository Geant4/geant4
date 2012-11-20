cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test16)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test16 test16.cc src/*.cc)
GEANT4_EXECUTABLE(test16_hadronic_exerciser test16_hadronic_exerciser.cc)
add_dependencies(test16 test16_hadronic_exerciser) # Artificial to force serialization (problem on windows)


#---Test definitions----------------------------------------------------------------------------------
GEANT4_ADD_TEST(test16-build BUILD test16)
GEANT4_ADD_TEST(test16_hadronic_exerciser COMMAND test16_hadronic_exerciser
                                          OUTPUT  test16.in
                                          DEPENDS test16-build)
GEANT4_ADD_TEST(test16 COMMAND test16 test16.in
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                       DEPENDS test16_hadronic_exerciser)
