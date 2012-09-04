GEANT4_ADD_UNIT_TESTS(testMatC*.cc testA*.cc testE*.cc testG*.cc
                      INCLUDE_DIRS geometry/management/include
                                   global/management/include
                                   global/HEPRandom/include
                                   global/HEPNumerics/include
                                   global/HEPGeometry/include
                                   materials/include
                                   intercoms/include
                                   ${CLHEP_INCLUDE_DIRS}
                      LIBRARIES    G4geometry)

