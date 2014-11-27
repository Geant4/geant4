#=============================================================================
# This file contains the definition of the tests with the GNUmake build
# It is only a sub-set of the existing examples to test that the build
# is correctly done. A couple of functions are included to facilititate
# the definition of the tests.
# Pere Mato 24/10/2012
#=============================================================================


#include(CMakeParseArguments)
#=============================================================================
# Copyright 2010 Alexander Neundorf <neundorf@kde.org>
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

function(CMAKE_PARSE_ARGUMENTS prefix _optionNames _singleArgNames _multiArgNames)
  # first set all result variables to empty/FALSE
  foreach(arg_name ${_singleArgNames} ${_multiArgNames})
    set(${prefix}_${arg_name})
  endforeach(arg_name)

  foreach(option ${_optionNames})
    set(${prefix}_${option} FALSE)
  endforeach(option)

  set(${prefix}_UNPARSED_ARGUMENTS)

  set(insideValues FALSE)
  set(currentArgName)

  # now iterate over all arguments and fill the result variables
  foreach(currentArg ${ARGN})
    list(FIND _optionNames "${currentArg}" optionIndex)  # ... then this marks the end of the arguments belonging to this keyword
    list(FIND _singleArgNames "${currentArg}" singleArgIndex)  # ... then this marks the end of the arguments belonging to this keyword
    list(FIND _multiArgNames "${currentArg}" multiArgIndex)  # ... then this marks the end of the arguments belonging to this keyword

    if(${optionIndex} EQUAL -1  AND  ${singleArgIndex} EQUAL -1  AND  ${multiArgIndex} EQUAL -1)
      if(insideValues)
        if("${insideValues}" STREQUAL "SINGLE")
          set(${prefix}_${currentArgName} ${currentArg})
          set(insideValues FALSE)
        elseif("${insideValues}" STREQUAL "MULTI")
          list(APPEND ${prefix}_${currentArgName} ${currentArg})
        endif()
      else(insideValues)
        list(APPEND ${prefix}_UNPARSED_ARGUMENTS ${currentArg})
      endif(insideValues)
    else()
      if(NOT ${optionIndex} EQUAL -1)
        set(${prefix}_${currentArg} TRUE)
        set(insideValues FALSE)
      elseif(NOT ${singleArgIndex} EQUAL -1)
        set(currentArgName ${currentArg})
        set(${prefix}_${currentArgName})
        set(insideValues "SINGLE")
      elseif(NOT ${multiArgIndex} EQUAL -1)
        set(currentArgName ${currentArg})
        set(${prefix}_${currentArgName})
        set(insideValues "MULTI")
      endif()
    endif()

  endforeach(currentArg)

  # propagate the result variables to the caller:
  foreach(arg_name ${_singleArgNames} ${_multiArgNames} ${_optionNames})
    set(${prefix}_${arg_name}  ${${prefix}_${arg_name}} PARENT_SCOPE)
  endforeach(arg_name)
  set(${prefix}_UNPARSED_ARGUMENTS ${${prefix}_UNPARSED_ARGUMENTS} PARENT_SCOPE)

endfunction(CMAKE_PARSE_ARGUMENTS _options _singleArgs _multiArgs)


#----------------------------------------------------------------------------
# function GEANT4_ADD_TEST( <name> COMMAND cmd [arg1... ] 
#                           [ENVIRONMENT var1=val1 var2=val2 ...
#                           [DEPENDS test1 ...]
#                           [TIMEOUT seconds] 
#                           [PASSREGEX exp] [FAILREGEX epx]
#                           [LABELS label1 label2 ...])
#
function(GEANT4_ADD_TEST test)
  CMAKE_PARSE_ARGUMENTS(ARG "" "TIMEOUT;PASSREGEX;FAILREGEX;BUILD;TARGET" "PRECMD;COMMAND;ENVIRONMENT;DEPENDS;LABELS" ${ARGN})
  if(ARG_BUILD)
    add_test(${test}-build make -C ${ARG_BUILD} ${ARG_TARGET})
    set(ARG_DEPENDS ${ARG_DEPENDS} ${test}-build)
    if(ARG_PRECMD)
        add_test(${test} ${ARG_PRECMD} && ${ARG_COMMAND})
    else()
        add_test(${test} ${ARG_COMMAND})
    endif()
  else()
    if(ARG_PRECMD)
        add_test(${test} ${ARG_PRECMD} && ${ARG_COMMAND})
    else()
        add_test(${test} ${ARG_COMMAND})
    endif()
  endif()
  if(ARG_ENVIRONMENT)
    #set_property(TEST ${test} PROPERTY ENVIRONMENT ${ARG_ENVIRONMENT})
    set(properties ${properties} ENVIRONMENT ${ARG_ENVIRONMENT})
  endif()
  if(ARG_TIMEOUT)
    set(properties ${properties} TIMEOUT ${ARG_TIMEOUT})
  endif()
  if(ARG_DEPENDS)
    set(properties ${properties} DEPENDS ${ARG_DEPENDS})
  endif()
  if(ARG_PASSREGEX)
    set(properties ${properties} PASS_REGULAR_EXPRESSION ${ARG_PASSREGEX})
  endif()
  if(ARG_FAILREGEX)
    set(properties ${properties} FAIL_REGULAR_EXPRESSION ${ARG_FAILREGEX})
   endif()
  if(ARG_LABELS)
    set(properties ${properties} LABELS ${ARG_LABELS})
  else()
    set(properties ${properties} LABELS Nightly)  
  endif()
  set_tests_properties(${test} PROPERTIES ${properties})
endfunction()


#-------------------------------------------------------------------------------

set(BINDIR $ENV{G4WORKDIR}/bin/$ENV{G4SYSTEM})
set(SRCDIR $ENV{G4INSTALL}/examples)

GEANT4_ADD_TEST(example-bas-b1  COMMAND ${BINDIR}/exampleB1 ${SRCDIR}/basic/B1/exampleB1.in
                                BUILD ${SRCDIR}/basic/B1)
GEANT4_ADD_TEST(example-bas-b2a COMMAND ${BINDIR}/exampleB2a ${SRCDIR}/basic/B2/B2a/exampleB2.in
                                BUILD ${SRCDIR}/basic/B2/B2a)
GEANT4_ADD_TEST(example-bas-b2b COMMAND ${BINDIR}/exampleB2b ${SRCDIR}/basic/B2/B2b/exampleB2.in
                                BUILD ${SRCDIR}/basic/B2/B2b)
GEANT4_ADD_TEST(example-bas-b3  COMMAND ${BINDIR}/exampleB3 ${SRCDIR}/basic/B3/exampleB3.in
                                BUILD ${SRCDIR}/basic/B3)
GEANT4_ADD_TEST(example-bas-b4a COMMAND ${BINDIR}/exampleB4a -m ${SRCDIR}/basic/B4/B4a/exampleB4.in
                                BUILD ${SRCDIR}/basic/B4/B4a)
GEANT4_ADD_TEST(example-bas-b4b COMMAND ${BINDIR}/exampleB4b -m ${SRCDIR}/basic/B4/B4b/exampleB4.in
                                BUILD ${SRCDIR}/basic/B4/B4b)
GEANT4_ADD_TEST(example-bas-b4c COMMAND ${BINDIR}/exampleB4c -m ${SRCDIR}/basic/B4/B4c/exampleB4.in
                                BUILD ${SRCDIR}/basic/B4/B4c)
GEANT4_ADD_TEST(example-bas-b4d COMMAND ${BINDIR}/exampleB4d -m ${SRCDIR}/basic/B4/B4d/exampleB4.in
                                BUILD ${SRCDIR}/basic/B4/B4d)
GEANT4_ADD_TEST(example-bas-b5  COMMAND ${BINDIR}/exampleB5 ${SRCDIR}/basic/B5/exampleB5.in
                                BUILD ${SRCDIR}/basic/B5)

foreach(_i 01 03) 
  GEANT4_ADD_TEST(example-ext-analysis-anaex${_i}-setup COMMAND make -C ${SRCDIR}/extended/analysis/AnaEx${_i} setup)
  GEANT4_ADD_TEST(example-ext-analysis-anaex${_i}-build COMMAND make -C ${SRCDIR}/extended/analysis/AnaEx${_i} DEPENDS example-ext-analysis-anaex${_i}-setup)
  GEANT4_ADD_TEST(example-ext-analysis-anaex${_i}-clean COMMAND make -C ${SRCDIR}/extended/analysis/AnaEx${_i} clean_setup DEPENDS example-ext-analysis-anaex${_i}-build)
  GEANT4_ADD_TEST(example-ext-analysis-anaex${_i} COMMAND ${BINDIR}/AnaEx${_i} ${SRCDIR}/extended/analysis/AnaEx${_i}/AnaEx${_i}.in DEPENDS example-ext-analysis-anaex${_i}-build)
endforeach()

GEANT4_ADD_TEST(example-ext-analysis-b1con COMMAND ${BINDIR}/exampleB1Con ${SRCDIR}/extended/analysis/B1Con/exampleB1Con.in
                                BUILD ${SRCDIR}/extended/analysis/B1Con)


foreach(_i 01 02 03) 
  GEANT4_ADD_TEST(example-ext-biasing-b${_i} COMMAND ${BINDIR}/exampleB${_i} 
                                        BUILD ${SRCDIR}/extended/biasing/B${_i})
endforeach()
foreach(_i 01 02 03 04) 
  GEANT4_ADD_TEST(example-ext-biasing-gb${_i} 
                COMMAND ${BINDIR}/exampleGB${_i} ${SRCDIR}/extended/biasing/GB${_i}/exampleGB${_i}.in
                BUILD ${SRCDIR}/extended/biasing/GB${_i})
endforeach()
GEANT4_ADD_TEST(example-ext-biasing-reversemc01 COMMAND ${BINDIR}/exampleRMC01 ${SRCDIR}/extended/biasing/ReverseMC01/run_adjoint_simulation_electron.g4mac  
                                        BUILD ${SRCDIR}/extended/biasing/ReverseMC01)
foreach(_i 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18)
  GEANT4_ADD_TEST(example-ext-electromagnetic-testem${_i}
                  COMMAND ${BINDIR}/TestEm${_i} ${SRCDIR}/extended/electromagnetic/TestEm${_i}/TestEm${_i}.in
                  BUILD ${SRCDIR}/extended/electromagnetic/TestEm${_i})
endforeach()

#GEANT4_ADD_TEST(example-ext-errorpropagation 
#                COMMAND ${BINDIR}/errprop
#                BUILD ${SRCDIR}/extended/errorpropagation
#        ------> ENVIRONMENT G4ERROR_TARGET=PLANE_SURFACE
#                            G4ERROR_MODE=BACKWARDS
#                            G4ERROR_PROP=UNTIL_TARGET)

GEANT4_ADD_TEST(example-ext-eventgenerator-exgps 
                COMMAND ${BINDIR}/exgps ${SRCDIR}/extended/eventgenerator/exgps/exgps_batch.in
                BUILD ${SRCDIR}/extended/eventgenerator/exgps)
GEANT4_ADD_TEST(example-ext-eventgenerator-particleGun-build 
                COMMAND make -C ${SRCDIR}/extended/eventgenerator/particleGun)
foreach(_i 1 2 3 4)
  GEANT4_ADD_TEST(example-ext-eventgenerator-particleGun-r${_i} 
                  COMMAND ${BINDIR}/particleGun ${SRCDIR}/extended/eventgenerator/particleGun/run${_i}.mac
                  DEPENDS example-ext-eventgenerator-particleGun-build)
endforeach()
GEANT4_ADD_TEST(example-ext-exoticphysics-monopole 
                COMMAND ${BINDIR}/monopole ${SRCDIR}/extended/exoticphysics/monopole/monopole.in
                BUILD ${SRCDIR}/extended/exoticphysics/monopole)

#GEANT4_ADD_TEST(example-ext-exoticphysics-phonons
#                COMMAND ${BINDIR}/XGeBox ${SRCDIR}/extended/exoticphysics/phonon/run.in
#                BUILD ${SRCDIR}/extended/exoticphysics/phonon
#       ------>  ENVIRONMENT CRYSTALMAPS=${SRCDIR}/extended/exoticphysics/phonon/CrystalMaps)

GEANT4_ADD_TEST(example-ext-exoticphysics-ucn
                 COMMAND ${BINDIR}/ExUCN -m ${SRCDIR}/extended/exoticphysics/ucn/ExUCN.in
                 BUILD ${SRCDIR}/extended/exoticphysics/ucn)

foreach(_i 01 02 03 04 05)
  GEANT4_ADD_TEST(example-ext-field-field${_i} 
                  COMMAND ${BINDIR}/field${_i} ${SRCDIR}/extended/field/field${_i}/field${_i}.in
                  BUILD ${SRCDIR}/extended/field/field${_i}) 
endforeach()
  GEANT4_ADD_TEST(example-ext-field-field06
                  COMMAND ${BINDIR}/field06 -m ${SRCDIR}/extended/field/field06/field06.in
                  BUILD ${SRCDIR}/extended/field/field06)
if(GEANT4_USE_G3TOG4)
  GEANT4_ADD_TEST(example-ext-g3tog4-clgeometry 
                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/g3tog4/clGeometry/clGeometry
                          ${SRCDIR}/extended/g3tog4/clGeometry/data/testmodel.dat
                          -m ${SRCDIR}/extended/g3tog4/clGeometry/clGeometry.in
                  BUILD ${SRCDIR}/extended/g3tog4/clGeometry
                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/g3tog4/clGeometry
                  BUILD clGeometry ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
endif()
foreach(_i 00 01 02 03 04 05 06)
  GEANT4_ADD_TEST(example-ext-hadronic-hadr${_i} 
                  COMMAND ${BINDIR}/Hadr${_i} ${SRCDIR}/extended/hadronic/Hadr${_i}/hadr${_i}.in
                  BUILD ${SRCDIR}/extended/hadronic/Hadr${_i})
endforeach()

#configure_file(${CTEST_SOURCE_DIRECTORY}/examples/extended/medical/DICOM/Data.dat ${CTEST_BINARY_DIRECTORY}/examples/extended/medical/DICOM/Data.dat)
#GEANT4_ADD_TEST(example-ext-medical-dicom 
#                COMMAND ${BINDIR}/DICOM ${SRCDIR}/extended/medical/DICOM/run.mac
#                BUILD ${SRCDIR}/extended/medical/DICOM)

#GEANT4_ADD_TEST(example-ext-medical-electronScattering 
#                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/medical/electronScattering/electronScattering
#                        ${SRCDIR}/extended/medical/electronScattering/electronScattering.in
#                PRECMD  ${CMAKE_COMMAND} -E make_directory results
#                BUILD ${SRCDIR}/extended/medical/electronScattering
#                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/medical/electronScattering
#                BUILD electronScattering ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-ext-medical-electronScattering2
                COMMAND ${BINDIR}/electronScattering2 ${SRCDIR}/extended/medical/electronScattering2/electronScattering2.in
                BUILD ${SRCDIR}/extended/medical/electronScattering2)
GEANT4_ADD_TEST(example-ext-medical-fanoCavity 
                COMMAND ${BINDIR}/fanoCavity ${SRCDIR}/extended/medical/fanoCavity/fanoCavity.in
                BUILD ${SRCDIR}/extended/medical/fanoCavity)
GEANT4_ADD_TEST(example-ext-medical-fanoCavity2 
                COMMAND ${BINDIR}/fanoCavity2
                        ${SRCDIR}/extended/medical/fanoCavity2/fanoCavity2.in
                BUILD ${SRCDIR}/extended/medical/fanoCavity2)
GEANT4_ADD_TEST(example-ext-medical-GammaTherapy 
                COMMAND ${BINDIR}/GammaTherapy ${SRCDIR}/extended/medical/GammaTherapy/GammaTherapy.in
                BUILD ${SRCDIR}/extended/medical/GammaTherapy)

set(DNA_EXAMPLE_SRC_DIR ${SRCDIR}/extended/medical/dna)

GEANT4_ADD_TEST(example-ext-medical-dna-dnaphysics 
                COMMAND ${BINDIR}/dnaphysics ${DNA_EXAMPLE_SRC_DIR}/dnaphysics/dnaphysics.in
                BUILD ${DNA_EXAMPLE_SRC_DIR}/dnaphysics)
 
GEANT4_ADD_TEST(example-ext-medical-dna-microdosimetry
                COMMAND ${BINDIR}/microdosimetry -mac ${DNA_EXAMPLE_SRC_DIR}/microdosimetry/microdosimetry.in
                BUILD ${DNA_EXAMPLE_SRC_DIR}/microdosimetry)

foreach(_i 1 2 3)                
# flag -mt always ON if G4MULTITHREADED is not defined, it will
# run in sequential mode 
GEANT4_ADD_TEST(example-ext-medical-dna-chem${_i}
                COMMAND ${BINDIR}/chem${_i} -mac ${DNA_EXAMPLE_SRC_DIR}/chem${_i}/beam.in
                BUILD ${DNA_EXAMPLE_SRC_DIR}/chem${_i})
endforeach()

#GEANT4_ADD_TEST(example-ext-medical-dna-wholeNuclearDNA
#                COMMAND ${BINDIR}/wholeNuclearDNA
#                        -mac ${DNA_EXAMPLE_SRC_DIR}/wholeNuclearDNA/wholeNuclearDNA.in
#                BUILD ${DNA_EXAMPLE_SRC_DIR}/wholeNuclearDNA)

#configure_file(${CTEST_SOURCE_DIRECTORY}examples/extended/medical/dna/pdb4dna/pdb4dna.in ${CTEST_BINARY_DIRECTORY}examples/extended/medical/dna/pdb4dna/pdb4dna.in)
#GEANT4_ADD_TEST(example-ext-medical-dna-pdb4dna
#                COMMAND ${BINDIR}/pdb4dna
#                        -mac ${DNA_EXAMPLE_SRC_DIR}/pdb4dna/pdb4dna.in
#                BUILD ${DNA_EXAMPLE_SRC_DIR}/pdb4dna)

GEANT4_ADD_TEST(example-ext-optical-opnovice 
                COMMAND ${BINDIR}/OpNovice -m ${SRCDIR}/extended/optical/OpNovice/OpNovice.in
                BUILD ${SRCDIR}/extended/optical/OpNovice)
GEANT4_ADD_TEST(example-ext-optical-lxe 
                COMMAND ${BINDIR}/LXe ${SRCDIR}/extended/optical/LXe/LXe.in
                BUILD ${SRCDIR}/extended/optical/LXe)
GEANT4_ADD_TEST(example-ext-optical-wls 
                COMMAND ${BINDIR}/wls ${SRCDIR}/extended/optical/wls/wls.in
                BUILD ${SRCDIR}/extended/optical/wls)
GEANT4_ADD_TEST(example-ext-parameterisations-gflash 
                COMMAND ${BINDIR}/ExGflash ${SRCDIR}/extended/parameterisations/gflash/test.mac
                BUILD ${SRCDIR}/extended/parameterisations/gflash)
GEANT4_ADD_TEST(example-ext-parameterisations-par01 COMMAND ${BINDIR}/examplePar01 ${SRCDIR}/extended/parameterisations/Par01/examplePar01.in
                BUILD ${SRCDIR}/extended/parameterisations/Par01)
#if(GEANT4_USE_GDML)
#  GEANT4_ADD_TEST(example-ext-persistency-gdml-g01 
#                  PRECMD  ${CMAKE_COMMAND} -E remove -f g01.gdml
#                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/persistency/gdml/G01/load_gdml
#                          ${SRCDIR}/extended/persistency/gdml/G01/solids.gdml 
#                          g01.gdml
#                          ${SRCDIR}/extended/persistency/gdml/G01/g01.in
#                  BUILD ${SRCDIR}/extended/persistency/gdml/G01
#                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/persistency/gdml/G01
#                  BUILD load_gdml)
#  GEANT4_ADD_TEST(example-ext-persistency-gdml-g04 
#                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/persistency/gdml/G04/gdml_det
#                          ${SRCDIR}/extended/persistency/gdml/G04/auxiliary.gdml 
#                          ${SRCDIR}/extended/persistency/gdml/G04/g04.mac
#                  BUILD ${SRCDIR}/extended/persistency/gdml/G04
#                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/persistency/gdml/G04
#                  BUILD gdml_det)
#endif()
#if(ROOT_FOUND AND BUILD_SHARED_LIBS)
#  GEANT4_ADD_TEST(example-ext-persistency-p01-write 
#                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/persistency/P01/exampleP01
#                          ${SRCDIR}/extended/persistency/P01/run.mac
#                  BUILD ${SRCDIR}/extended/persistency/P01
#                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/persistency/P01
#                  BUILD exampleP01 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
#                  
#  GEANT4_ADD_TEST(example-ext-persistency-p01-read 
#                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/persistency/P01/readHits hits.root
#                  BUILD ${SRCDIR}/extended/persistency/P01
#                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/persistency/P01
#                  BUILD readHits ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT}
#                  DEPENDS example-ext-persistency-p01-write)
#endif()

#GEANT4_ADD_TEST(example-ext-persistency-p03 
#                PRECMD  ${CMAKE_COMMAND} -E copy ${SRCDIR}/extended/persistency/P03/g4geom.txt g4geom.txt
#                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/persistency/P03/textGeom
#                        ${SRCDIR}/extended/persistency/P03/batch.mac
#                BUILD ${SRCDIR}/extended/persistency/P03
#                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/persistency/P03
#                BUILD textGeom ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
#GEANT4_ADD_TEST(example-ext-polarisation-pol01 
#                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/polarisation/Pol01/pol01
#                        ${SRCDIR}/extended/polarisation/Pol01/pol01.in
#                BUILD ${SRCDIR}/extended/polarisation/Pol01
#                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/polarisation/Pol01
#                BUILD pol01)
#foreach(_i 01 02)
#  GEANT4_ADD_TEST(example-ext-radioactivedecay-rdecay${_i} 
#                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/radioactivedecay/rdecay${_i}/rdecay${_i}
#                        ${SRCDIR}/extended/radioactivedecay/rdecay${_i}/rdecay${_i}.in
#                BUILD ${SRCDIR}/extended/radioactivedecay/rdecay${_i}
#                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/radioactivedecay/rdecay${_i}
#                BUILD rdecay${_i}  ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
#endforeach()
#GEANT4_ADD_TEST(example-ext-runandevent-re01 
#                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE01/RE01
#                        ${SRCDIR}/extended/runAndEvent/RE01/sample.in
#                ERROR sample.err
#                PRECMD  ${CMAKE_COMMAND} -E copy ${SRCDIR}/extended/runAndEvent/RE01/pythia_event.data pythia_event.data
#                POSTCMD ${CMAKE_COMMAND} -E compare_files sample.err ${SRCDIR}/extended/runAndEvent/RE01/sample.err
#                BUILD ${SRCDIR}/extended/runAndEvent/RE01
#                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE01
#                BUILD RE01 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
#GEANT4_ADD_TEST(example-ext-runandevent-re02 
#                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02/exampleRE02
#                        ${SRCDIR}/extended/runAndEvent/RE02/run.mac
#                BUILD ${SRCDIR}/extended/runAndEvent/RE02
#                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02
#                BUILD exampleRE02 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
#GEANT4_ADD_TEST(example-ext-runandevent-re02-run3 
#                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02/exampleRE02
#                        ${SRCDIR}/extended/runAndEvent/RE02/run3.mac
#                BUILD ${SRCDIR}/extended/runAndEvent/RE02
#                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02
#                DEPENDS example-ext-runandevent-re02 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
#GEANT4_ADD_TEST(example-ext-runandevent-re02-run4 
#                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02/exampleRE02
#                        ${SRCDIR}/extended/runAndEvent/RE02/run4.mac
#                BUILD ${SRCDIR}/extended/runAndEvent/RE02
#                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE02
#                DEPENDS example-ext-runandevent-re02 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
#file(WRITE ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE04/score.ref "G4VScoringMesh::DrawColorChart(): no visualization system\n")  
#GEANT4_ADD_TEST(example-ext-runandevent-re04 
#                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE04/RE04
#                        ${SRCDIR}/extended/runAndEvent/RE04/score.mac
#                ERROR score.err
#                POSTCMD ${CMAKE_COMMAND} -E compare_files score.err score.ref
#                BUILD ${SRCDIR}/extended/runAndEvent/RE04
#                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/runAndEvent/RE04
#                BUILD RE04 ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
GEANT4_ADD_TEST(example-ext-runandevent-re05
                COMMAND ${BINDIR}/exampleRE05 ${SRCDIR}/extended/runAndEvent/RE05/exampleRE05.in
                BUILD ${SRCDIR}/extended/runAndEvent/RE05)
GEANT4_ADD_TEST(example-ext-runandevent-re06 
                COMMAND ${BINDIR}/exampleRE06 ${SRCDIR}/extended/runAndEvent/RE06/exampleRE06.in
                BUILD ${SRCDIR}/extended/runAndEvent/RE06)
#GEANT4_ADD_TEST(example-ext-visualization-uservisaction-build
#                BUILD ${SRCDIR}/extended/visualization/userVisAction
#                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/visualization/userVisAction
#                BUILD userVisAction)
#GEANT4_ADD_TEST(example-ext-visualization-uservisaction-uva
#                COMMAND ${CMAKE_BINARY_DIR}/examples/extended/visualization/userVisAction/userVisAction
#                        ${SRCDIR}/extended/visualization/userVisAction/userVisAction.in
#                SOURCE_DIR ${SRCDIR}/extended/visualization/userVisAction
#                BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/visualization/userVisAction
#                DEPENDS example-ext-visualization-uservisaction-build ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
#foreach(_i 1 2)
#  GEANT4_ADD_TEST(example-ext-visualization-uservisaction-run${_i}
#                  COMMAND ${CMAKE_BINARY_DIR}/examples/extended/visualization/userVisAction/userVisAction
#                          ${SRCDIR}/extended/visualization/userVisAction/run${_i}.mac
#                  SOURCE_DIR ${SRCDIR}/extended/visualization/userVisAction
#                  BINARY_DIR ${CMAKE_BINARY_DIR}/examples/extended/visualization/userVisAction
#                  DEPENDS example-ext-visualization-uservisaction-build ENVIRONMENT ${GEANT4_TEST_ENVIRONMENT})
#endforeach()
#

#-------------------------------------------------------------------------------------------------


