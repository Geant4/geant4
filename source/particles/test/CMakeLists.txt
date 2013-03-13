
configure_file(checkParticles/mass_width_2011.mcd mass_width_2011.mcd COPYONLY)

GEANT4_ADD_UNIT_TESTS(checkParticles printMassTable reportParticleTable
                      INCLUDE_DIRS global/management/include
                                   global/HEPRandom/include
                                   global/HEPNumerics/include
                                   global/HEPGeometry/include
                                   track/include
                                   materials/include
                                   intercoms/include
                                   processes/management/include
                                   particles/management/include
                                   particles/bosons/include
                                   particles/leptons/include
                                   particles/hadrons/mesons/include 
                                   particles/hadrons/barions/include
                                   particles/hadrons/ions/include
                                   particles/shortlived/include
                                   particles/utils/include
                                   ${CLHEP_INCLUDE_DIRS}
                      LIBRARIES    G4particles)

