cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

#---For IDEs the location of executables is prepended with the configuration (Release, Debug,...)
if(CMAKE_GENERATOR MATCHES Makefiles)
  set(cfg)
else()
  set(cfg $<CONFIGURATION>/)
endif()

#---Need to see if ROOT is available to enable some of the tests
find_package(ROOT QUIET)
find_package(AIDA QUIET)
find_package(gdml QUIET)

#-------------------------------------------------------------------------------------------------
#---Consistency-----------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
#GEANT4_ADD_TEST(validate-sources 
#                COMMAND ${CMAKE_MAKE_PROGRAM} validate_sources
#                WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                LABELS Continuous) 

#-------------------------------------------------------------------------------------------------
#---Basic-----------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
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
                        -m ${CMAKE_SOURCE_DIR}/examples/basic/B4/B4a/exampleB4.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/basic/B4/B4a
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/basic/B4/B4a
                BUILD exampleB4a PROJECT B4a ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-bas-b4b 
                COMMAND ${CMAKE_BINARY_DIR}/examples/basic/B4/B4b/exampleB4b
                        -m ${CMAKE_SOURCE_DIR}/examples/basic/B4/B4b/exampleB4.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/basic/B4/B4b
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/basic/B4/B4b
                BUILD exampleB4b PROJECT B4b ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-bas-b4c 
                COMMAND ${CMAKE_BINARY_DIR}/examples/basic/B4/B4c/exampleB4c
                        -m ${CMAKE_SOURCE_DIR}/examples/basic/B4/B4c/exampleB4.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/basic/B4/B4c
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/basic/B4/B4c
                BUILD exampleB4c PROJECT B4c ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-bas-b4d 
                COMMAND ${CMAKE_BINARY_DIR}/examples/basic/B4/B4d/exampleB4d
                        -m ${CMAKE_SOURCE_DIR}/examples/basic/B4/B4d/exampleB4.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/basic/B4/B4d
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/basic/B4/B4d
                BUILD exampleB4d PROJECT B4d ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-bas-b5 
                COMMAND ${CMAKE_BINARY_DIR}/examples/basic/B5/exampleB5
                        ${CMAKE_SOURCE_DIR}/examples/basic/B5/exampleB5.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/basic/B5
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/basic/B5
                BUILD exampleB5 PROJECT B5 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

#-------------------------------------------------------------------------------------------------
#---Extended-----------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
GEANT4_ADD_TEST(example-ext-analysis-anaex01 
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/analysis/AnaEx01/AnaEx01 
                        ${CMAKE_SOURCE_DIR}/examples/extended/analysis/AnaEx01/AnaEx01.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/analysis/AnaEx01
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/analysis/AnaEx01
                BUILD AnaEx01 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
if(ROOT_FOUND)
  GEANT4_ADD_TEST(example-ext-analysis-anaex02 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/analysis/AnaEx02/AnaEx02 
                          ${CMAKE_SOURCE_DIR}/examples/extended/analysis/AnaEx02/AnaEx02.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/analysis/AnaEx02
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/analysis/AnaEx02
                  BUILD AnaEx02 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
if(AIDA_FOUND)
  GEANT4_ADD_TEST(example-ext-analysis-anaex03
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/analysis/AnaEx03/AnaEx03
                          ${CMAKE_SOURCE_DIR}/examples/extended/analysis/AnaEx03/AnaEx03.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/analysis/AnaEx03
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/analysis/AnaEx03
                  BUILD AnaEx03 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
GEANT4_ADD_TEST(example-ext-analysis-b1con 
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/analysis/B1Con/exampleB1Con 
                        ${CMAKE_SOURCE_DIR}/examples/extended/analysis/B1Con/exampleB1Con.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/analysis/B1Con
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/analysis/B1Con
                BUILD exampleB1Con PROJECT B1Con ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
foreach(_i 01 02 03)
  GEANT4_ADD_TEST(example-ext-biasing-b${_i} 
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/biasing/B${_i}/exampleB${_i}
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/biasing/B${_i}
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/biasing/B${_i}
                BUILD exampleB${_i} PROJECT B${_i} ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()
foreach(_i 01 02 03 04)
  GEANT4_ADD_TEST(example-ext-biasing-gb${_i} 
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/biasing/GB${_i}/exampleGB${_i}
                        ${CMAKE_SOURCE_DIR}/examples/extended/biasing/GB${_i}/exampleGB${_i}.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/biasing/GB${_i}
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/biasing/GB${_i}
                BUILD exampleGB${_i} PROJECT GB${_i}project ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()
GEANT4_ADD_TEST(example-ext-biasing-reversemc01 
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/biasing/ReverseMC01/exampleRMC01
                        ${CMAKE_SOURCE_DIR}/examples/extended/biasing/ReverseMC01/run_adjoint_simulation_electron.mac
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/biasing/ReverseMC01
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/biasing/ReverseMC01
                BUILD exampleRMC01 PROJECT ReverseMC01 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
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
                BUILD errprop PROJECT errorpropagation 
                ENVIRONMENT G4ERROR_TARGET=PLANE_SURFACE
                            G4ERROR_MODE=BACKWARDS
                            G4ERROR_PROP=UNTIL_TARGET )
GEANT4_ADD_TEST(example-ext-eventgenerator-exgps 
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/eventgenerator/exgps/exgps
                        ${CMAKE_SOURCE_DIR}/examples/extended/eventgenerator/exgps/exgps_batch.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/eventgenerator/exgps
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/eventgenerator/exgps
                BUILD exgps 
		ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
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
                BUILD monopole ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
		)
include( ${CMAKE_SOURCE_DIR}/examples/extended/exoticphysics/phonon/CTestDefinitions.txt)
foreach(_i 01 02 03 04 05)
  GEANT4_ADD_TEST(example-ext-field-field${_i} 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/field/field${_i}/field${_i} 
                          ${CMAKE_SOURCE_DIR}/examples/extended/field/field${_i}/field${_i}.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/field/field${_i} 
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/field/field${_i} 
                  BUILD field${_i} ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()
  GEANT4_ADD_TEST(example-ext-field-field06
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/field/field06/field06
                          -m ${CMAKE_SOURCE_DIR}/examples/extended/field/field06/field06.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/field/field06
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/field/field06
                  BUILD field06 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
if(GEANT4_USE_G3TOG4)
  GEANT4_ADD_TEST(example-ext-g3tog4-clgeometry 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/g3tog4/clGeometry/clGeometry
                          ${CMAKE_SOURCE_DIR}/examples/extended/g3tog4/clGeometry/data/testmodel.dat
                          -m ${CMAKE_SOURCE_DIR}/examples/extended/g3tog4/clGeometry/clGeometry.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/g3tog4/clGeometry
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/g3tog4/clGeometry
                  BUILD clGeometry ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
GEANT4_ADD_TEST(example-ext-geometry-transforms 
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/geometry/transforms/transforms
                        ${CMAKE_SOURCE_DIR}/examples/extended/geometry/transforms/debug.mac
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/geometry/transforms
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/geometry/transforms
                BUILD transforms )
foreach(_i 00 01 02 03 04 05 06)
  GEANT4_ADD_TEST(example-ext-hadronic-hadr${_i} 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/hadronic/Hadr${_i}/Hadr${_i}
                          ${CMAKE_SOURCE_DIR}/examples/extended/hadronic/Hadr${_i}/hadr${_i}.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/hadronic/Hadr${_i}
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/hadronic/Hadr${_i}
                  BUILD Hadr${_i} ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()
GEANT4_ADD_TEST(example-ext-hadronic-FissionFragment 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/hadronic/FissionFragment/FissionFragment
						        -i ${CMAKE_SOURCE_DIR}/examples/extended/hadronic/FissionFragment/batch.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/hadronic/FissionFragment
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/hadronic/FissionFragment
                  BUILD FissionFragment 
						ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} 
						         G4NEUTRON_HP_USE_WENDT_FISSION_MODEL=1
							      G4NEUTRONHP_PRODUCE_FISSION_FRAGMENTS=1)
GEANT4_ADD_TEST(example-ext-medical-dicom 
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/medical/DICOM/DICOM
                        ${CMAKE_SOURCE_DIR}/examples/extended/medical/DICOM/run.mac
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/medical/DICOM
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/medical/DICOM
                BUILD DICOM ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-ext-medical-electronScattering 
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/medical/electronScattering/electronScattering
                        ${CMAKE_SOURCE_DIR}/examples/extended/medical/electronScattering/electronScattering.in
                PRECMD  ${CMAKE_COMMAND} -E make_directory results
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/medical/electronScattering
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/medical/electronScattering
                BUILD electronScattering ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-ext-medical-electronScattering2 
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/medical/electronScattering2/electronScattering2
                        ${CMAKE_SOURCE_DIR}/examples/extended/medical/electronScattering2/electronScattering2.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/medical/electronScattering2
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/medical/electronScattering2
                BUILD electronScattering2 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
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
GEANT4_ADD_TEST(example-ext-optical-opnovice
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/optical/OpNovice/OpNovice 
                        -m ${CMAKE_SOURCE_DIR}/examples/extended/optical/OpNovice/OpNovice.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/optical/OpNovice
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/optical/OpNovice
                BUILD OpNovice ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
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
                BUILD ExGflash PROJECT gflash ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-ext-parameterisations-par01 
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/parameterisations/Par01/examplePar01 
                        ${CMAKE_SOURCE_DIR}/examples/extended/parameterisations/Par01/examplePar01.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/parameterisations/Par01
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/parameterisations/Par01
                BUILD examplePar01 PROJECT Par01 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
if(GEANT4_USE_GDML)
  GEANT4_ADD_TEST(example-ext-persistency-gdml-g01 
                  PRECMD  ${CMAKE_COMMAND} -E remove -f g01.gdml
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/persistency/gdml/G01/load_gdml
                          ${CMAKE_SOURCE_DIR}/examples/extended/persistency/gdml/G01/solids.gdml 
                          g01.gdml
                          ${CMAKE_SOURCE_DIR}/examples/extended/persistency/gdml/G01/g01.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/persistency/gdml/G01
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/persistency/gdml/G01
                  BUILD load_gdml PROJECT G01 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
  GEANT4_ADD_TEST(example-ext-persistency-gdml-g04 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/persistency/gdml/G04/gdml_det
                          ${CMAKE_SOURCE_DIR}/examples/extended/persistency/gdml/G04/auxiliary.gdml 
                          ${CMAKE_SOURCE_DIR}/examples/extended/persistency/gdml/G04/g04.mac
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/persistency/gdml/G04
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/persistency/gdml/G04
                  BUILD gdml_det PROJECT G04 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
if(ROOT_FOUND AND GCCXML AND Geant4_shared_FOUND)
  GEANT4_ADD_TEST(example-ext-persistency-p01-write 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/persistency/P01/exampleP01
                          ${CMAKE_SOURCE_DIR}/examples/extended/persistency/P01/run.mac
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/persistency/P01
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/persistency/P01
                  BUILD exampleP01 PROJECT P01 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
                  
  GEANT4_ADD_TEST(example-ext-persistency-p01-read 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/persistency/P01/readHits hits.root
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/persistency/P01
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/persistency/P01
                  BUILD readHits PROJECT P01 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                  DEPENDS example-ext-persistency-p01-write)
endif()

GEANT4_ADD_TEST(example-ext-persistency-p03 
                PRECMD  ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/examples/extended/persistency/P03/g4geom.txt g4geom.txt
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/persistency/P03/textGeom
                        ${CMAKE_SOURCE_DIR}/examples/extended/persistency/P03/batch.mac
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/persistency/P03
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/persistency/P03
                BUILD textGeom  PROJECT P03 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-ext-polarisation-pol01 
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/polarisation/Pol01/Pol01
                        ${CMAKE_SOURCE_DIR}/examples/extended/polarisation/Pol01/pol01.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/polarisation/Pol01
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/polarisation/Pol01
                BUILD Pol01 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
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
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02/RE02
                        ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE02/run.mac
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE02
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02
                BUILD RE02 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-ext-runandevent-re02-run3 
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02/RE02
                        ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE02/run3.mac
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE02
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02
                DEPENDS example-ext-runandevent-re02 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-ext-runandevent-re02-run4 
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02/RE02
                        ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE02/run4.mac
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE02
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02
                DEPENDS example-ext-runandevent-re02 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
file(WRITE ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE04/score.ref "G4VScoringMesh::DrawColorChart(): no visualization system\n")  
GEANT4_ADD_TEST(example-ext-runandevent-re04 
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE04/RE04
                        ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE04/score.mac
                ERROR score.err
                POSTCMD ${CMAKE_COMMAND} -E compare_files score.err score.ref
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE04
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE04
                BUILD RE04 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
                LABELS Nightly Continuous)
GEANT4_ADD_TEST(example-ext-runandevent-re05
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE05/exampleRE05
                        ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE05/exampleRE05.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE05
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE05
                BUILD exampleRE05 PROJECT RE05 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-ext-runandevent-re06
                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE06/exampleRE06
                        ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE06/exampleRE06.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/runAndEvent/RE06
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE06
                BUILD exampleRE06 PROJECT RE06 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
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
include( ${CMAKE_SOURCE_DIR}/examples/extended/exoticphysics/channeling/CTestDefinitions.txt)

GEANT4_ADD_TEST(example-ext-exoticphysics-ExUCN
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/exoticphysics/ucn/ExUCN
                          -m ${CMAKE_SOURCE_DIR}/examples/extended/exoticphysics/ucn/ExUCN.in
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/exoticphysics/ucn
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/exoticphysics/ucn
                  BUILD ExUCN ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

set(DNA_EXAMPLE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/extended/medical/dna)
set(DNA_EXAMPLE_BIN_DIR ${CMAKE_BINARY_DIR}/examples/extended/medical/dna)

GEANT4_ADD_TEST(example-ext-medical-dna-pdb4dna
                COMMAND ${DNA_EXAMPLE_BIN_DIR}/pdb4dna/pdb4dna
                        -mac ${DNA_EXAMPLE_SOURCE_DIR}/pdb4dna/pdb4dna.in
                SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/pdb4dna
                BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/pdb4dna
                BUILD pdb4dna ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

GEANT4_ADD_TEST(example-ext-medical-dna-wholeNuclearDNA
                COMMAND ${DNA_EXAMPLE_BIN_DIR}/wholeNuclearDNA/wholeNuclearDNA 
                        -mac ${DNA_EXAMPLE_SOURCE_DIR}/wholeNuclearDNA/wholeNuclearDNA.in
                SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/wholeNuclearDNA
                BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/wholeNuclearDNA
                BUILD wholeNuclearDNA ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
                 
GEANT4_ADD_TEST(example-ext-medical-dna-dnaphysics 
                COMMAND ${DNA_EXAMPLE_BIN_DIR}/dnaphysics/dnaphysics 
                        ${DNA_EXAMPLE_SOURCE_DIR}/dnaphysics/dnaphysics.in
                SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/dnaphysics
                BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/dnaphysics
                BUILD dnaphysics ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
                 
GEANT4_ADD_TEST(example-ext-medical-dna-microdosimetry 
                COMMAND ${DNA_EXAMPLE_BIN_DIR}/microdosimetry/microdosimetry 
                        -mac ${DNA_EXAMPLE_SOURCE_DIR}/microdosimetry/microdosimetry.in
                SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/microdosimetry
                BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/microdosimetry
                BUILD microdosimetry ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
                
foreach(_i 1 2 3)          
# flag -mt always ON if G4MULTITHREADED is not defined, it will
# run in sequential mode      
GEANT4_ADD_TEST(example-ext-medical-dna-chem${_i}
                COMMAND ${DNA_EXAMPLE_BIN_DIR}/chem${_i}/chem${_i} 
                        -mac ${DNA_EXAMPLE_SOURCE_DIR}/chem${_i}/beam.in 
                SOURCE_DIR ${DNA_EXAMPLE_SOURCE_DIR}/chem${_i}
                BINARY_DIR ${DNA_EXAMPLE_BIN_DIR}/chem${_i}
                BUILD chem${_i} ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endforeach()

#-------------------------------------------------------------------------------------------------
#---Advanced--------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------
GEANT4_ADD_TEST(example-adv-air_shower 
                COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/air_shower/Ultra
                        ${CMAKE_SOURCE_DIR}/examples/advanced/air_shower/UltraMacro.mac
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/air_shower
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/air_shower
                BUILD Ultra PROJECT air_shower
		ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-adv-amsEcal 
                COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/amsEcal/amsEcal 
                        ${CMAKE_SOURCE_DIR}/examples/advanced/amsEcal/run1.mac
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/amsEcal
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/amsEcal
                BUILD amsEcal ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-adv-brachytherapy 
                COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/brachytherapy/Brachy 
                        ${CMAKE_SOURCE_DIR}/examples/advanced/brachytherapy/IridiumSourceMacro.mac
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/brachytherapy
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/brachytherapy
                BUILD Brachy PROJECT Brachy
		ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
# ChargeExchangeMC has a mandatory GDML dependency
if (gdml_FOUND)
 GEANT4_ADD_TEST(example-adv-ChargeExchangeMC-build
                 SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/ChargeExchangeMC
                 BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/ChargeExchangeMC
                 BUILD ChargeExchangeMC PROJECT ChargeExchangeMC)
endif() 
GEANT4_ADD_TEST(example-adv-composite_calorimeter-build
                 SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/composite_calorimeter
                 BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/composite_calorimeter
                 BUILD CompositeCalorimeter PROJECT CompositeCalorimeter)
GEANT4_ADD_TEST(example-adv-eRosita 
                COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/eRosita/eRosita 
                        ${CMAKE_SOURCE_DIR}/examples/advanced/eRosita/eRosita.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/eRosita
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/eRosita
                BUILD eRosita ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-adv-gammaknife
                COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/gammaknife/gammaknife 
                        ${CMAKE_SOURCE_DIR}/examples/advanced/gammaknife/simulation.mac
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/gammaknife
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/gammaknife
                BUILD gammaknife  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-adv-gammaray_telescope 
                COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/gammaray_telescope/GammaRayTel 
                        ${CMAKE_SOURCE_DIR}/examples/advanced/gammaray_telescope/gammaraytel.in
#                ERROR gammaraytel.err
#                POSTCMD ${CMAKE_COMMAND} -E compare_files gammaraytel.err ${CMAKE_SOURCE_DIR}/examples/advanced/gammaray_telescope/gammaraytel.err
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/gammaray_telescope
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/gammaray_telescope
                BUILD GammaRayTel  PROJECT gammaray_telescope
		ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-adv-microbeam
                COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/microbeam/microbeam
                        ${CMAKE_SOURCE_DIR}/examples/advanced/microbeam/microbeam.mac
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/microbeam
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/microbeam
                BUILD microbeam ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-adv-microelectronics 
                COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/microelectronics/microelectronics 
                        ${CMAKE_SOURCE_DIR}/examples/advanced/microelectronics/microelectronics.in
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/microelectronics
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/microelectronics
                BUILD microelectronics ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
if(GEANT4_USE_SYSTEM_CLHEP)
  GEANT4_ADD_TEST(example-adv-nanobeam 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/nanobeam/nanobeam 
                          ${CMAKE_SOURCE_DIR}/examples/advanced/nanobeam/default.mac
                  SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/nanobeam
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/nanobeam
                  BUILD nanobeam  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
GEANT4_ADD_TEST(example-adv-hadrontherapy 
                COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/hadrontherapy/hadrontherapy
                        ${CMAKE_SOURCE_DIR}/examples/advanced/hadrontherapy/batch.mac
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/hadrontherapy
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/hadrontherapy
                BUILD hadrontherapy ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-adv-human_phantom 
                COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/human_phantom/phantom batch.mac
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/human_phantom
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/human_phantom
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/examples/advanced/human_phantom
                BUILD phantom PROJECT human_phantom ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-adv-purging_magnet
                COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/purging_magnet/PurgMag MacroTesting.mac
                WORKING_DIRECTORY  ${CMAKE_SOURCE_DIR}/examples/advanced/purging_magnet
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/purging_magnet
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/purging_magnet
                BUILD PurgMag PROJECT purging_magnet ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-adv-radioprotection
                COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/radioprotection/radioprotection
                        ${CMAKE_SOURCE_DIR}/examples/advanced/radioprotection/run.mac
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/radioprotection
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/radioprotection
                BUILD radioprotection ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-adv-xray_fluorescence
                COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/xray_fluorescence/XrayFluo livermore.mac 1
                WORKING_DIRECTORY  ${CMAKE_SOURCE_DIR}/examples/advanced/xray_fluorescence
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/xray_fluorescence
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/xray_fluorescence
                BUILD XrayFluo PROJECT xray_fluorescence ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT} XRAYDATA=. )
GEANT4_ADD_TEST(example-adv-xray_telescope
                COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/xray_telescope/XrayTel test.mac
                WORKING_DIRECTORY  ${CMAKE_SOURCE_DIR}/examples/advanced/xray_telescope
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/xray_telescope
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/xray_telescope
                BUILD XrayTel PROJECT xray_telescope ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-adv-iort_therapy-build
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/iort_therapy
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/iort_therapy
                BUILD iort_therapy)
GEANT4_ADD_TEST(example-adv-medical_linac-build
                SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/medical_linac
                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/medical_linac
                BUILD ml2 PROJECT medical_linac)
GEANT4_ADD_TEST(example-adv-underground_physics
		COMMAND ${CMAKE_BINARY_DIR}/examples/advanced/underground_physics/DMX
			${CMAKE_BINARY_DIR}/examples/advanced/underground_physics/gamma_1000.mac
		SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/underground_physics
		BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced/underground_physics
		BUILD DMX PROJECT underground_physics
		ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-adv-lArCal
		COMMAND
                ${CMAKE_BINARY_DIR}/examples/advanced/lAr_calorimeter/lArCal
			${CMAKE_BINARY_DIR}/examples/advanced/lAr_calorimeter/run1.mac
		SOURCE_DIR ${CMAKE_SOURCE_DIR}/examples/advanced/lAr_calorimeter
		BINARY_DIR ${CMAKE_BINARY_DIR}/examples/advanced//lAr_calorimeter
		BUILD lArCal 
		ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})

#-------------------------------------------------------------------------------------------------
#---Simplified Calorimeter and other PhysicsChecks tests----------------------------------------------------
#-------------------------------------------------------------------------------------------------
include( ${CMAKE_SOURCE_DIR}/verification/SimplifiedCalorimeter/physval/CTestDefinitions.txt OPTIONAL)
include( ${CMAKE_SOURCE_DIR}/examples/extended/eventgenerator/exgps/CTestsDefinition.txt OPTIONAL)

#-------------------------------------------------------------------------------------------------
#---CTests for MIC (Intel Xeon Phi)
#-------------------------------------------------------------------------------------------------
include( ${CMAKE_SOURCE_DIR}/tests/ctests/CTests-mic.txt OPTIONAL)

