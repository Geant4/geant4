
GEANT4_ADD_UNIT_TESTS(INCLUDE_DIRS geometry/management/include
                                   geometry/solids/CSG/include
                                   geometry/solids/specific/include
                                   global/management/include
                                   global/HEPRandom/include
                                   global/HEPNumerics/include
                                   global/HEPGeometry/include
                                   graphics_reps/include
                                   ${CLHEP_INCLUDE_DIRS}
                      LIBRARIES G4geometry)
  


