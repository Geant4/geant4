# Copy selected reference outputs from test to source area ready to be
#  committed from any suitable system (currently sungeant with export G4DEBUG=1).
#  Edit g4_src below to select your own checked out directory.
#   Run ". setup.sh" first to set variables G4INSTALL, G4SYSTEM and G4WORKDIR.
#   Then
#
#  cvs commit -m "New output" tests examples
#
# Run . setup.sh first to set variables G4INSTALL, G4SYSTEM and G4WORKDIR
#
#g4_test=$G4WORKDIR/stt.geant4-05-00-ref-03+tags2038/$G4SYSTEM
g4_test=$G4WORKDIR/stt/$G4SYSTEM
g4_src=/afs/cern.ch/rd44/user/mclareni/geant4
#
cp $g4_test/test01.out           $g4_src/tests/test01/test01.out
cp $g4_test/test02.hadron.out    $g4_src/tests/test02/test02.hadron.out
cp $g4_test/test02.out           $g4_src/tests/test02/test02.out
cp $g4_test/test05.out           $g4_src/tests/test05/test05.out
cp $g4_test/test07.out           $g4_src/tests/test07/test07.out
cp $g4_test/test09.out           $g4_src/tests/test09/test09.out
cp $g4_test/test10.out           $g4_src/tests/test10/test10.out
cp $g4_test/test11.out           $g4_src/tests/test11/test11.out
cp $g4_test/test12.out           $g4_src/tests/test12/test12.out
cp $g4_test/test13.out           $g4_src/tests/test13/test13.out
cp $g4_test/test14.out           $g4_src/tests/test14/test14.out
cp $g4_test/test15.out           $g4_src/tests/test15/test15.out
cp $g4_test/test16.out           $g4_src/tests/test16/test16.out
cp $g4_test/test17.out           $g4_src/tests/test17/test17.out
cp $g4_test/test18.out           $g4_src/tests/test18/test18.out
cp $g4_test/test19.out           $g4_src/tests/test19/test19.out
cp $g4_test/test20.out           $g4_src/tests/test20/test20.out
cp $g4_test/test21.out           $g4_src/tests/test21/test21.out
cp $g4_test/test22.out           $g4_src/tests/test22/test22.out
cp $g4_test/test23.out           $g4_src/tests/test23/test23.out
cp $g4_test/test24.out           $g4_src/tests/test24/test24.out
cp $g4_test/test25.out           $g4_src/tests/test25/test25.out
# cp $g4_test/test26.out           $g4_src/tests/test26/test26.out
cp $g4_test/test27.out           $g4_src/tests/test27/test27.out
cp $g4_test/test28.out           $g4_src/tests/test28/test28.out
cp $g4_test/test29.out           $g4_src/tests/test29/test29.out
cp $g4_test/test32.out           $g4_src/tests/test32/test32.out
cp $g4_test/test33.out           $g4_src/tests/test33/MassGeo_TimedApp.out
cp $g4_test/test33_1.out         $g4_src/tests/test33/ParallelGeo_TimedApp.out
cp $g4_test/test34.out           $g4_src/tests/test34/test34.out
cp $g4_test/test39.out           $g4_src/tests/test39/test39.out
cp $g4_test/test40.out           $g4_src/tests/test40/test40.out

cp $g4_test/test101.out          $g4_src/examples/novice/N01/exampleN01.out
cp $g4_test/test102.out          $g4_src/examples/novice/N02/exampleN02.out
cp $g4_test/test103.out          $g4_src/examples/novice/N03/exampleN03.out
cp $g4_test/test104.EMtest.out   $g4_src/examples/novice/N04/exampleN04.EMtest.out
cp $g4_test/test104.out          $g4_src/examples/novice/N04/exampleN04.out
cp $g4_test/test105.out          $g4_src/examples/novice/N05/exampleN05.out
cp $g4_test/test106.out          $g4_src/examples/novice/N06/exampleN06.out  
cp $g4_test/test107.out          $g4_src/examples/novice/N07/exampleN07.out

cp $g4_test/test500.out          $g4_src/examples/extended/electromagnetic/TestEm0/TestEm0.out
cp $g4_test/test501.out          $g4_src/examples/extended/electromagnetic/TestEm1/TestEm1.out
cp $g4_test/test5010.out         $g4_src/examples/extended/electromagnetic/TestEm10/TestEm10.out
cp $g4_test/test5011.out         $g4_src/examples/extended/electromagnetic/TestEm11/TestEm11.out
cp $g4_test/test5012.out         $g4_src/examples/extended/electromagnetic/TestEm12/TestEm12.out
cp $g4_test/test5013.out         $g4_src/examples/extended/electromagnetic/TestEm13/TestEm13.out
cp $g4_test/test5014.out         $g4_src/examples/extended/electromagnetic/TestEm14/TestEm14.out
cp $g4_test/test5015.out         $g4_src/examples/extended/electromagnetic/TestEm15/TestEm15.out
cp $g4_test/test5016.out         $g4_src/examples/extended/electromagnetic/TestEm16/TestEm16.out
cp $g4_test/test5017.out         $g4_src/examples/extended/electromagnetic/TestEm17/TestEm17.out
cp $g4_test/test5018.out         $g4_src/examples/extended/electromagnetic/TestEm18/TestEm18.out
cp $g4_test/test502.out          $g4_src/examples/extended/electromagnetic/TestEm2/TestEm2.out
cp $g4_test/test503.out          $g4_src/examples/extended/electromagnetic/TestEm3/TestEm3.out
cp $g4_test/test504.out          $g4_src/examples/extended/electromagnetic/TestEm4/TestEm4.out
cp $g4_test/test505.out          $g4_src/examples/extended/electromagnetic/TestEm5/TestEm5.out  
cp $g4_test/test506.out          $g4_src/examples/extended/electromagnetic/TestEm6/TestEm6.out
cp $g4_test/test507.out          $g4_src/examples/extended/electromagnetic/TestEm7/TestEm7.out
cp $g4_test/test508.out          $g4_src/examples/extended/electromagnetic/TestEm8/TestEm8.out
cp $g4_test/test509.out          $g4_src/examples/extended/electromagnetic/TestEm9/TestEm9.out

cp $g4_test/test601.out          $g4_src/examples/extended/g3tog4/clGeometry/clGeometry.out
cp $g4_test/test602.out          $g4_src/examples/extended/g3tog4/cltog4/cltog4.out
cp $g4_test/test701.out          $g4_src/examples/extended/field/field01/field01.out
cp $g4_test/test702.out          $g4_src/examples/extended/field/field02/field02.out
cp $g4_test/test703.out          $g4_src/examples/extended/field/field03/field03.out
cp $g4_test/test704.out          $g4_src/examples/extended/field/field04/field04.out
cp $g4_test/test801.out          $g4_src/examples/extended/biasing/B01/exampleB01.out
cp $g4_test/test1002.out         $g4_src/examples/advanced/microbeam/microbeam.out
cp $g4_test/test1003.out         $g4_src/examples/advanced/raredecay_calorimetry/PhotIn.out
cp $g4_test/test1004.out         $g4_src/examples/advanced/nanobeam/nanobeam.out
cp $g4_test/test1008.out         $g4_src/examples/advanced/cosmicray_charging/shoot.out
cp $g4_test/test2001.out         $g4_src/examples/extended/medical/GammaTherapy/GammaTherapy.out
cp $g4_test/test2002.out         $g4_src/examples/extended/medical/fanoCavity/fanoCavity.out
cp $g4_test/test2003.out         $g4_src/examples/extended/medical/fanoCavity2/fanoCavity2.out
cp $g4_test/test6001.out         $g4_src/examples/extended/hadronic/Hadr01/hadr01.out
cp $g4_test/test7001.out         $g4_src/examples/extended/polarisation/Pol01/pol01.out
