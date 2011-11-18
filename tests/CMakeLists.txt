#---Add all subdirectories that contains a CMakeLists.txt file -------------------------------------
file(GLOB_RECURSE files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} CMakeLists.txt)
foreach( file ${files} )
  get_filename_component(path ${file} PATH)
  if(path)
    add_subdirectory(${path})
  endif()
endforeach()


#---Add tests for the examples-----------------------------------------------------------------------
#---Basic-----------------------------------------------------------------------------------------
  GEANT4_ADD_TEST(example-bas-b1 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/basic/B1/exampleB1 
                          ${CMAKE_SOURCE_DIR}/examples/basic/B1/exampleB1.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/basic/B1
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/basic/B1
                  BUILD exampleB1 PROJECT B1 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-bas-b2a 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/basic/B2/B2a/exampleB2a
                          ${CMAKE_SOURCE_DIR}/examples/basic/B2/B2a/exampleB2.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/basic/B2/B2a
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/basic/B2/B2a
                  BUILD exampleB2a PROJECT B2a ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-bas-b2b 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/basic/B2/B2b/exampleB2b
                          ${CMAKE_SOURCE_DIR}/examples/basic/B2/B2b/exampleB2.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/basic/B2/B2b
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/basic/B2/B2b
                  BUILD exampleB2b PROJECT B2b ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-bas-b3 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/basic/B3/exampleB3 
                          ${CMAKE_SOURCE_DIR}/examples/basic/B3/exampleB3.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/basic/B3
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/basic/B3
                  BUILD exampleB3 PROJECT B3 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-bas-b4a 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/basic/B4/B4a/exampleB4a
                          ${CMAKE_SOURCE_DIR}/examples/basic/B4/B4a/exampleB4.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/basic/B4/B4a
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/basic/B4/B4a
                  BUILD exampleB4a PROJECT B4a ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-bas-b4b 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/basic/B4/B4b/exampleB4b
                          ${CMAKE_SOURCE_DIR}/examples/basic/B4/B4b/exampleB4.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/basic/B4/B4b
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/basic/B4/B4b
                  BUILD exampleB4b PROJECT B4b ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-bas-b4c 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/basic/B4/B4c/exampleB4c
                          ${CMAKE_SOURCE_DIR}/examples/basic/B4/B4c/exampleB4.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/basic/B4/B4c
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/basic/B4/B4c
                  BUILD exampleB4c PROJECT B4c ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-bas-b4d 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/basic/B4/B4d/exampleB4d
                          ${CMAKE_SOURCE_DIR}/examples/basic/B4/B4d/exampleB4.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/basic/B4/B4d
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/basic/B4/B4d
                  BUILD exampleB4d PROJECT B4d ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

  #---Novice-----------------------------------------------------------------------------------------
  GEANT4_ADD_TEST(example-nov-n01 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/novice/N01/exampleN01 
                          ${CMAKE_SOURCE_DIR}/examples/novice/N01/exampleN01.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/novice/N01
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/novice/N01
                  BUILD exampleN01 PROJECT N01)
  GEANT4_ADD_TEST(example-nov-n02 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/novice/N02/exampleN02
                          ${CMAKE_SOURCE_DIR}/examples/novice/N02/exampleN02.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/novice/N02
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/novice/N02
                  BUILD exampleN02 PROJECT N02 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-nov-n03 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/novice/N03/exampleN03 
                          ${CMAKE_SOURCE_DIR}/examples/novice/N03/exampleN03.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/novice/N03
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/novice/N03
                  BUILD exampleN03 PROJECT N03 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-nov-n04 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/novice/N04/exampleN04 
                          ${CMAKE_SOURCE_DIR}/examples/novice/N04/exampleN04.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/novice/N04
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/novice/N04
                  BUILD exampleN04 PROJECT N04
                  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-nov-n05 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/novice/N05/exampleN05 
                          ${CMAKE_SOURCE_DIR}/examples/novice/N05/exampleN05.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/novice/N05
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/novice/N05
                  BUILD exampleN05 PROJECT N05 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-nov-n06 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/novice/N06/exampleN06 
                          ${CMAKE_SOURCE_DIR}/examples/novice/N06/exampleN06.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/novice/N06
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/novice/N06
                  BUILD exampleN06 PROJECT N06 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-nov-n07 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/novice/N07/exampleN07 
                          ${CMAKE_SOURCE_DIR}/examples/novice/N07/exampleN07.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/novice/N07
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/novice/N07
                  BUILD exampleN07 PROJECT N07 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  #---Extended-----------------------------------------------------------------------------------------
  GEANT4_ADD_TEST(example-ext-analysis-a01 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/analysis/A01/A01app 
                          ${CMAKE_SOURCE_DIR}/examples/extended/analysis/A01/test.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/analysis/A01
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/analysis/A01
                  BUILD A01app PROJECT A01 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-ext-analysis-anaex01 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/analysis/AnaEx01/AnaEx01 
                          ${CMAKE_SOURCE_DIR}/examples/extended/analysis/AnaEx01/AnaEx01.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/analysis/AnaEx01
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/analysis/AnaEx01
                  BUILD AnaEx01 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-ext-analysis-anaex02 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/analysis/AnaEx02/AnaEx02 
                          ${CMAKE_SOURCE_DIR}/examples/extended/analysis/AnaEx02/AnaEx02.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/analysis/AnaEx02
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/analysis/AnaEx02
                  BUILD AnaEx02 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-ext-analysis-n03con 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/analysis/N03Con/exampleN03Con 
                          ${CMAKE_SOURCE_DIR}/examples/extended/analysis/N03Con/exampleN03Con.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/analysis/N03Con
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/analysis/N03Con
                  BUILD exampleN03Con ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-ext-biasing-b01 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/biasing/B01/exampleB01
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/biasing/B01
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/biasing/B01
                  BUILD exampleB01 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-ext-biasing-b02 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/biasing/B02/exampleB02
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/biasing/B02
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/biasing/B02
                  BUILD exampleB02 PROJECT B02)
  GEANT4_ADD_TEST(example-ext-biasing-reversemc01 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/biasing/ReverseMC01/exampleRMC01
                          ${CMAKE_SOURCE_DIR}/examples/extended/biasing/ReverseMC01/run_adjoint_simulation_electron.g4mac
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/biasing/ReverseMC01
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/biasing/ReverseMC01
                  BUILD exampleRMC01 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  foreach(_i 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18)
    GEANT4_ADD_TEST(example-ext-electromagnetic-testem${_i}
                    COMMAND ${CMAKE_BINARY_DIR}/examples/extended/electromagnetic/TestEm${_i}/TestEm${_i} 
                            ${CMAKE_SOURCE_DIR}/examples/extended/electromagnetic/TestEm${_i}/TestEm${_i}.in
                    SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/electromagnetic/TestEm${_i}
                    BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/electromagnetic/TestEm${_i}
                    BUILD TestEm${_i}
                    ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  endforeach()
  GEANT4_ADD_TEST(example-ext-errorpropagation 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/errorpropagation/errprop
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/errorpropagation
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/errorpropagation
                  BUILD errprop
                  ENVIRONMENT G4ERROR_TARGET=PLANE_SURFACE
                              G4ERROR_MODE=BACKWARDS
                              G4ERROR_PROP=UNTIL_TARGET )
  GEANT4_ADD_TEST(example-ext-eventgenerator-exgps 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/eventgenerator/exgps/exGPS
                          ${CMAKE_SOURCE_DIR}/examples/extended/eventgenerator/exgps/exgps_batch.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/eventgenerator/exgps
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/eventgenerator/exgps
                  BUILD exGPS)
  GEANT4_ADD_TEST(example-ext-eventgenerator-particleGun-build 
                    SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/eventgenerator/particleGun
                    BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/eventgenerator/particleGun
                    BUILD particleGun)
  foreach(_i 1 2 3 4)
    GEANT4_ADD_TEST(example-ext-eventgenerator-particleGun-r${_i} 
                    COMMAND ${CMAKE_BINARY_DIR}/examples/extended/eventgenerator/particleGun/particleGun 
                            ${CMAKE_SOURCE_DIR}/examples/extended/eventgenerator/particleGun/run${_i}.mac
                    SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/eventgenerator/particleGun
                    BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/eventgenerator/particleGun
                    DEPENDS example-ext-eventgenerator-particleGun-build)
  endforeach()
  GEANT4_ADD_TEST(example-ext-exoticphysics-monopole 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/exoticphysics/monopole/monopole
                          ${CMAKE_SOURCE_DIR}/examples/extended/exoticphysics/monopole/monopole.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/exoticphysics/monopole
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/exoticphysics/monopole
                  BUILD monopole ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} )
  foreach(_i 01 02 03 04 05)
    GEANT4_ADD_TEST(example-ext-field-field${_i} 
                    COMMAND ${CMAKE_BINARY_DIR}/examples/extended/field/field${_i}/field${_i} 
                            ${CMAKE_SOURCE_DIR}/examples/extended/field/field${_i}/field${_i}.in
                    SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/field/field${_i} 
                    BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/field/field${_i} 
                    BUILD field${_i} ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  endforeach()
  if(GEANT4_USE_G3TOG4)
    GEANT4_ADD_TEST(example-ext-g3tog4-clgeometry 
                    COMMAND ${CMAKE_BINARY_DIR}/examples/extended/g3tog4/clGeometry/clGeometry
                            ${CMAKE_SOURCE_DIR}/examples/extended/g3tog4/data/testmodel.dat
                            ${CMAKE_SOURCE_DIR}/examples/extended/g3tog4/clGeometry/clGeometry.in
                    SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/g3tog4/clGeometry
                    BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/g3tog4/clGeometry
                    BUILD clGeometry ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
    GEANT4_ADD_TEST(example-ext-g3tog4-cltog4 
                    COMMAND ${CMAKE_BINARY_DIR}/examples/extended/g3tog4/cltog4/cltog4
                            ${CMAKE_SOURCE_DIR}/examples/extended/g3tog4/data/testmodel.dat
                            ${CMAKE_SOURCE_DIR}/examples/extended/g3tog4/cltog4/cltog4.in
                    SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/g3tog4/cltog4
                    BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/g3tog4/cltog4
                    BUILD cltog4 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  endif()
  GEANT4_ADD_TEST(example-ext-geometry-olap 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/geometry/olap/olapex
                          ${CMAKE_SOURCE_DIR}/examples/extended/geometry/olap/batch.mac
                  ERROR batch.err
                  POSTCMD ${CMAKE_COMMAND} -E compare_files batch.err ${CMAKE_SOURCE_DIR}/examples/extended/geometry/olap/batch.err
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/geometry/olap
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/geometry/olap
                  BUILD olapex)
  foreach(_i 00 01 02)
    GEANT4_ADD_TEST(example-ext-hadronic-hadr${_i} 
                    COMMAND ${CMAKE_BINARY_DIR}/examples/extended/hadronic/Hadr${_i}/hadr${_i}
                            ${CMAKE_SOURCE_DIR}/examples/extended/hadronic/Hadr${_i}/hadr${_i}.in
                    SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/hadronic/Hadr${_i}
                    BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/hadronic/Hadr${_i}
                    BUILD hadr${_i} ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  endforeach()
  GEANT4_ADD_TEST(example-ext-medical-dicom 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/medical/DICOM/dicom
                          ${CMAKE_SOURCE_DIR}/examples/extended/medical/DICOM/run.mac
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/medical/DICOM
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/medical/DICOM
                  BUILD dicom)
  GEANT4_ADD_TEST(example-ext-medical-electronScattering 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/medical/electronScattering/electronScattering
                          ${CMAKE_SOURCE_DIR}/examples/extended/medical/electronScattering/electronScattering.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/medical/electronScattering
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/medical/electronScattering
                  BUILD electronScattering ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-ext-medical-fanoCavity 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/medical/fanoCavity/fanoCavity
                          ${CMAKE_SOURCE_DIR}/examples/extended/medical/fanoCavity/fanoCavity.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/medical/fanoCavity
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/medical/fanoCavity
                  BUILD fanoCavity)
  GEANT4_ADD_TEST(example-ext-medical-fanoCavity2 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/medical/fanoCavity2/fanoCavity2
                          ${CMAKE_SOURCE_DIR}/examples/extended/medical/fanoCavity2/fanoCavity2.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/medical/fanoCavity2
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/medical/fanoCavity2
                  BUILD fanoCavity2)
  GEANT4_ADD_TEST(example-ext-medical-GammaTherapy 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/medical/GammaTherapy/GammaTherapy
                          ${CMAKE_SOURCE_DIR}/examples/extended/medical/GammaTherapy/GammaTherapy.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/medical/GammaTherapy
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/medical/GammaTherapy
                  BUILD GammaTherapy ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-ext-optical-lxe 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/optical/LXe/LXe
                          ${CMAKE_SOURCE_DIR}/examples/extended/optical/LXe/LXe.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/optical/LXe
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/optical/LXe
                  BUILD LXe ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-ext-optical-wls 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/optical/wls/wls
                          ${CMAKE_SOURCE_DIR}/examples/extended/optical/wls/wls.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/optical/wls
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/optical/wls
                  BUILD wls  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-ext-parameterisations-gflash 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/parameterisations/gflash/ExGflash
                          ${CMAKE_SOURCE_DIR}/examples/extended/parameterisations/gflash/test.mac
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/parameterisations/gflash
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/parameterisations/gflash
                  BUILD ExGflash ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  if(GEANT4_USE_GDML)
    GEANT4_ADD_TEST(example-ext-persistency-gdml-g01 
                    PRECMD  ${CMAKE_COMMAND} -E remove -f g01.gdml
                    COMMAND ${CMAKE_BINARY_DIR}/examples/extended/persistency/gdml/G01/load_gdml
                            ${CMAKE_SOURCE_DIR}/examples/extended/persistency/gdml/G01/solids.gdml 
                            g01.gdml
                            ${CMAKE_SOURCE_DIR}/examples/extended/persistency/gdml/G01/g01.in
                    SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/persistency/gdml/G01
                    BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/persistency/gdml/G01
                    BUILD load_gdml)
    GEANT4_ADD_TEST(example-ext-persistency-gdml-g04 
                    COMMAND ${CMAKE_BINARY_DIR}/examples/extended/persistency/gdml/G04/gdml_det
                            ${CMAKE_SOURCE_DIR}/examples/extended/persistency/gdml/G04/auxiliary.gdml 
                            ${CMAKE_SOURCE_DIR}/examples/extended/persistency/gdml/G04/g04.mac
                    SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/persistency/gdml/G04
                    BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/persistency/gdml/G04
                    BUILD gdml_det)
  endif()
  GEANT4_ADD_TEST(example-ext-persistency-p03 
                  PRECMD  ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/examples/extended/persistency/P03/g4geom.txt g4geom.txt
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/persistency/P03/textGeom
                          ${CMAKE_SOURCE_DIR}/examples/extended/persistency/P03/batch.mac
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/persistency/P03
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/persistency/P03
                  BUILD textGeom ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-ext-polarisation-pol01 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/polarisation/Pol01/pol01
                          ${CMAKE_SOURCE_DIR}/examples/extended/polarisation/Pol01/pol01.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/polarisation/Pol01
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/polarisation/Pol01
                  BUILD pol01)
  foreach(_i 01 02)
    GEANT4_ADD_TEST(example-ext-radioactivedecay-rdecay${_i} 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/radioactivedecay/rdecay${_i}/rdecay${_i}
                          ${CMAKE_SOURCE_DIR}/examples/extended/radioactivedecay/rdecay${_i}/rdecay${_i}.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/radioactivedecay/rdecay${_i}
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/radioactivedecay/rdecay${_i}
                  BUILD rdecay${_i}  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  endforeach()
  GEANT4_ADD_TEST(example-ext-runandevent-re01 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE01/RE01
                          ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE01/sample.in
                  ERROR sample.err
                  PRECMD  ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE01/pythia_event.data pythia_event.data
                  POSTCMD ${CMAKE_COMMAND} -E compare_files sample.err ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE01/sample.err
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE01
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE01
                  BUILD RE01 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-ext-runandevent-re02 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02/exampleRE02
                          ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE02/run.mac
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE02
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02
                  BUILD exampleRE02 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-ext-runandevent-re02-run3 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02/exampleRE02
                          ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE02/run3.mac
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE02
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02
                  DEPENDS example-ext-runandevent-re02 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-ext-runandevent-re02-run4 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02/exampleRE02
                          ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE02/run4.mac
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE02
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02
                  DEPENDS example-ext-runandevent-re02 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-ext-visualization-uservisaction-build
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/visualization/userVisAction
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/visualization/userVisAction
                  BUILD userVisAction)
  GEANT4_ADD_TEST(example-ext-visualization-uservisaction-uva
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/visualization/userVisAction/userVisAction
                          ${CMAKE_SOURCE_DIR}/examples/extended/visualization/userVisAction/userVisAction.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/visualization/userVisAction
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/visualization/userVisAction
                  DEPENDS example-ext-visualization-uservisaction-build ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  foreach(_i 1 2)
    GEANT4_ADD_TEST(example-ext-visualization-uservisaction-run${_i}
                    COMMAND ${CMAKE_BINARY_DIR}/examples/extended/visualization/userVisAction/userVisAction
                            ${CMAKE_SOURCE_DIR}/examples/extended/visualization/userVisAction/run${_i}.mac
                    SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/visualization/userVisAction
                    BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/visualization/userVisAction
                    DEPENDS example-ext-visualization-uservisaction-build ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  endforeach()


  #---Advanced-----------------------------------------------------------------------------------------------
  GEANT4_ADD_TEST(example-adv-amsEcal 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/amsEcal/amsEcal 
                          ${CMAKE_SOURCE_DIR}/examples/advanced/amsEcal/run1.mac
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/amsEcal
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/amsEcal
                  BUILD amsEcal ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-adv-dnaphysics 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/dnaphysics/DNAPhysics 
                          ${CMAKE_SOURCE_DIR}/examples/advanced/dnaphysics/dna.mac
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/dnaphysics
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/dnaphysics
                  BUILD DNAPhysics ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-adv-eRosita 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/eRosita/eRosita 
                          ${CMAKE_SOURCE_DIR}/examples/advanced/eRosita/eRosita.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/eRosita
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/eRosita
                  BUILD eRosita ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-adv-gammaray_telescope 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/gammaray_telescope/GammaRayTel 
                          ${CMAKE_SOURCE_DIR}/examples/advanced/gammaray_telescope/gammaraytel.in
                  ERROR gammaraytel.err
                  POSTCMD ${CMAKE_COMMAND} -E compare_files gammaraytel.err ${CMAKE_SOURCE_DIR}/examples/advanced/gammaray_telescope/gammaraytel.err
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/gammaray_telescope
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/gammaray_telescope
                  BUILD GammaRayTel  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-adv-microdosimetry 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/microdosimetry/Microdosimetry 
                          ${CMAKE_SOURCE_DIR}/examples/advanced/microdosimetry/microdosimetry.mac
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/microdosimetry
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/microdosimetry
                  BUILD Microdosimetry ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  if(GEANT4_USE_SYSTEM_CLHEP)
    GEANT4_ADD_TEST(example-adv-nanobeam 
                    COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/nanobeam/Nanobeam 
                            ${CMAKE_SOURCE_DIR}/examples/advanced/nanobeam/default.mac
                    SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/nanobeam
                    BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/nanobeam
                    BUILD Nanobeam  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  endif()


