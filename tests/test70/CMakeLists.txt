cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(test70)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

GEANT4_EXECUTABLE(test70 test70.cc src/*.cc)

#---Test definitions----------------------------------------------------------------------------------
configure_file(commands.mac commands.mac COPYONLY) # test expects file to be local
GEANT4_ADD_TEST(test70 COMMAND test70
                       BUILD test70
                       ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

