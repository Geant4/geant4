#----------------------------------------------------------------------------
# intercoms unit test programs
#
if(CMAKE_PROJECT_NAME STREQUAL Geant4)

  GEANT4_ADD_UNIT_TESTS(INCLUDE_DIRS intercoms/include
                                     global/management/include
                                     global/HEPRandom/include
                                     global/HEPNumerics/include
                                     global/HEPGeometry/include
                                     ${CLHEP_INCLUDE_DIRS}
                        LIBRARIES G4intercoms)
else()

  cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
  project(testGenericMessenger)
  find_package(Geant4 REQUIRED)
  include(${Geant4_USE_FILE})
  GEANT4_EXECUTABLE(testGenericMessenger testGenericMessenger.cc)

endif()
