GEANT4_ADD_UNIT_TESTS(test*.cc field07 
                      INCLUDE_DIRS geometry/management/include
                                   geometry/volumes/include
                                   geometry/navigation/include
                                   geometry/solids/CSG/include
                                   geometry/solids/Boolean/include
                                   geometry/solids/specific/include
                                   geometry/magneticfield/include
                                   global/management/include
                                   global/HEPRandom/include
                                   global/HEPNumerics/include
                                   global/HEPGeometry/include
                                   graphics_reps/include
                                   run/include
                                   particles/management/include
                                   particles/bosons/include
                                   particles/leptons/include
                                   particles/hadrons/mesons/include
                                   particles/hadrons/barions/include
                                   particles/hadrons/ions/include
                                   materials/include
                                   event/include
                                   digits_hits/hits/include
                                   digits_hits/digits/include
                                   digits_hits/detector/include
                                   tracking/include
                                   track/include
                                   processes/management/include
                                   processes/cuts/include
                                   processes/electromagnetic/standard/include
                                   processes/electromagnetic/utils/include
                                   processes/electromagnetic/muons/include
                                   processes/electromagnetic/highenergy/include
                                   processes/decay/include
                                   processes/transportation/include
                                   intercoms/include
                                   interfaces/basic/include
                                   interfaces/common/include
                                   ${CLHEP_INCLUDE_DIRS}
                      LIBRARIES G4geometry G4particles G4interfaces G4run)

