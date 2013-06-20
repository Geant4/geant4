#!/bin/sh
# Copy selected reference outputs from test to source area ready to be
#  committed from any suitable system.
#  Edit g4_src below to select your own checked out directory.
#   Run ". setup.sh" first to set variables G4INSTALL, G4SYSTEM and G4WORKDIR.
#   Then
#
#  use tagrefout script to commit and tag changed files

# to add new outputs:
# in bash do:
# for f in `svn st | grep ? | awk '{ print $2 }'`; do svn add $f;done
#
# tag :
tag=geant4-09-06-ref-05
svn_branch="`echo $tag | sed -e 's/geant4-//'`_branch"

# Patch/release  branch...
#svn_branch=geant4-09-06-patches_branch
#tag=geant4-09-06-patch-02_cand-01
#tag="`echo $svn_branch | sed -e 's/patches_branch/patch-02/'`"

slot=nightly
slot=release
platform=slc5-gcc43-RelWithDebInfo

#  use outputs from build with tag.... or ref tag see below
#g4_test=/build/cdash/....
#g4_test=/build/cdash/G4/$slot/$platform/

g4_test=/ec/G4-builds/$slot/$platform/build

#using refout in ~stesting:
g4_src=~/refout/$tag/geant4



echo "Start to check out directories test and examples"
echo "hit return to continue"
read ok

function update() {
	new=$g4_test/$1/$2.out
	#echo "		From= $new"
	dest=$g4_src/$1/$3
	#echo "		To  = $dest"
	if test -r $new ; then 
	   echo "            copy $new $dest"  
	   cp $new $dest
	else  
           echo "File $new not found"
	fi  
	#read bla 
}

if ( test -d $g4_src ) then
  echo "Error: geant4 source already exists"
  exit 1
fi  
mkdir -p $g4_src || exit

# was used with cvs
#cd $g4_src/..
#rmdir geant4
#cvs -d :gserver:stesting@geant4.cvs.cern.ch:/cvs/Geant4 co \
#    -r $tag -l geant4 > cvs.log 2>&1
# cd geant4
# cvs update -d -P tests examples > ../cvs.log 2>&1

cd $g4_src
svn co svn+ssh://svn.cern.ch/reps/geant4/branches/geant4/_symbols/$svn_branch/tests tests
svn co svn+ssh://svn.cern.ch/reps/geant4/branches/geant4/_symbols/$svn_branch/examples examples
 
echo "will copy reference ouputs from" 
echo "  $g4_test"
echo " to"
echo "  $g4_src"
echo " please check date of one test01 logfile:"
ls  -l $g4_test/tests/test01/test01.out
echo "hit return to continue"
read ok


update  tests/test01  test01      test01.out
update  tests/test02  test02-had  test02.hadron.out
update  tests/test02  test02      test02.out
update  tests/test05  test05      test05.out
update  tests/test07  test07      test07.out
update  tests/test09  test09      test09.out
update  tests/test10  test10      test10.out
update  tests/test11  test11      test11.out
update  tests/test12  test12-QGSP_FTFP_BERT      test12-QGSP_FTFP_BERT.out
update  tests/test12  test12      test12.out
update  tests/test13  test13      test13.out
update  tests/test14  test14-std      	test14_std.out
update  tests/test14  test14-lowE 		test14_lowE.out
update  tests/test14  test14-penelope	test14_penelope.out
update  tests/test15  test15      test15.out
update  tests/test16  test16      test16.out
update  tests/test17  test17      test17.out
update  tests/test18  test18      test18.out
update  tests/test19  test19      test19.out
# retired 9.3r? update  tests/test20  test20      test20.out
# retired 9.5r02 update  tests/test21  test21      test21.out
update  tests/test22  test22      test22.out
update  tests/test23  test23      test23.out
update  tests/test24  test24      test24.out
update  tests/test25  test25      test25.out
update  tests/test27  test27      test27.out
update  tests/test28  test28      test28.out
update  tests/test29  test29      test29.out
update  tests/test31  test31      test31.out
#retired 9.3.r09 update  tests/test32  test32      test32.out
update  tests/test33  test33      MassGeo_TimedApp.out
update  tests/test33  test33-1    ParallelGeo_TimedApp.out
update  tests/test34  test34      test34.out
update  tests/test35  test35      test35.out
update  tests/test37  test37      test37.out
update  tests/test39  test39      test39.out
update  tests/test40  test40      test40.out
update  tests/test41  test41      test41.out
update  tests/test44  test44      test44.out
update  tests/test49  test49      test49.out
update  tests/test53  test53      test53.out
update  tests/test54  test54      test54.out
update  tests/test55  test55      test55.out
update  tests/test58  test58      test58.out
update  tests/test60  test60      test60.out
update  tests/test61  test61      test61.out
update  tests/test62  test62      test62.out
update  tests/test67  test67      test67.out
update  tests/test70  test70      test70.out
update  tests/test74  test74      test74.out


update  examples/extended/analysis/A01     	example-ext-analysis-*	   test.out
update  examples/extended/analysis/N03Con     	example-ext-analysis-*	   exampleN03Con.out

update  examples/extended/biasing/B01/     	example-ext-biasing	   exampleB01.out
update  examples/extended/biasing/ReverseMC01/  example-ext-biasing	   run_adjoint_simulation_electron.out

update  examples/extended/electromagnetic/TestEm0      example-ext-electromagnetic-*       TestEm0.out
update  examples/extended/electromagnetic/TestEm1      example-ext-electromagnetic-*       TestEm1.out
update  examples/extended/electromagnetic/TestEm2      example-ext-electromagnetic-*       TestEm2.out
update  examples/extended/electromagnetic/TestEm3      example-ext-electromagnetic-*       TestEm3.out
update  examples/extended/electromagnetic/TestEm4      example-ext-electromagnetic-*       TestEm4.out
update  examples/extended/electromagnetic/TestEm5      example-ext-electromagnetic-*       TestEm5.out  
update  examples/extended/electromagnetic/TestEm6      example-ext-electromagnetic-*       TestEm6.out
update  examples/extended/electromagnetic/TestEm7      example-ext-electromagnetic-*       TestEm7.out
update  examples/extended/electromagnetic/TestEm8      example-ext-electromagnetic-*       TestEm8.out
update  examples/extended/electromagnetic/TestEm9      example-ext-electromagnetic-*       TestEm9.out
update  examples/extended/electromagnetic/TestEm10     example-ext-electromagnetic-*       TestEm10.out
update  examples/extended/electromagnetic/TestEm11     example-ext-electromagnetic-*       TestEm11.out
update  examples/extended/electromagnetic/TestEm12     example-ext-electromagnetic-*       TestEm12.out
update  examples/extended/electromagnetic/TestEm13     example-ext-electromagnetic-*       TestEm13.out
update  examples/extended/electromagnetic/TestEm14     example-ext-electromagnetic-*       TestEm14.out
update  examples/extended/electromagnetic/TestEm15     example-ext-electromagnetic-*       TestEm15.out
update  examples/extended/electromagnetic/TestEm16     example-ext-electromagnetic-*       TestEm16.out
update  examples/extended/electromagnetic/TestEm17     example-ext-electromagnetic-*       TestEm17.out
update  examples/extended/electromagnetic/TestEm18     example-ext-electromagnetic-*       TestEm18.out

update  examples/extended/errorpropagation		example-ext-errorpropagation	errorprop.out

update  examples/extended/eventgenerator/exgps		example-ext-errorpropagation   exgps_batch.out
#update  examples/extended/eventgenerator/particleGun	*-r1-run        run1.out
#update  examples/extended/eventgenerator/particleGun	*-r2-run        run2.out
#update  examples/extended/eventgenerator/particleGun	*-r3-run        run3.out
#update  examples/extended/eventgenerator/particleGun	*-r4-run        run4.out

update  examples/extended/exoticphysics/monopole	example-ext-exoticphysics-monopole	monopole.out

update  examples/extended/field/field01      		example-ext-field-*	 field01.out
update  examples/extended/field/field02      		example-ext-field-*	 field02.out
update  examples/extended/field/field03      		example-ext-field-*	 field03.out
update  examples/extended/field/field04      		example-ext-field-*	 field04.out
update  examples/extended/field/field05      		example-ext-field-*	 field05.out

#update  examples/extended/g3tog4/clGeometry          	*-run	   clGeometry.out
#update  examples/extended/g3tog4/cltog4      		*-run	   cltog4.out

update  examples/extended/geometry/olap     		example-ext-geometry-olap	   batch.out

update  examples/extended/hadronic/Hadr00	      	example-ext-hadronic-*	   hadr00.out
update  examples/extended/hadronic/Hadr01	      	example-ext-hadronic-*	   hadr01.out
update  examples/extended/hadronic/Hadr02	      	example-ext-hadronic-*	   hadr02.out
update  examples/extended/hadronic/Hadr03	      	example-ext-hadronic-*	   hadr03.out

update  examples/extended/medical/DICOM			example-ext-medical-*      run.out
update  examples/extended/medical/electronScattering    example-ext-medical-*      electronScattering.out
update  examples/extended/medical/fanoCavity      	example-ext-medical-*	   fanoCavity.out
update  examples/extended/medical/fanoCavity2      	example-ext-medical-*	   fanoCavity2.out
update  examples/extended/medical/GammaTherapy      	example-ext-medical-*	   GammaTherapy.out

update  examples/extended/optical/LXe			example-ext-optical-*	   LXe.out		
update  examples/extended/optical/wls			example-ext-optical-*	   wls.out

update  examples/extended/parameterisations/gflash     	example-ext-parameterisations-*	   test.out

update  examples/extended/polarisation/Pol01      	example-ext-polarisation-*	   pol01.out

update  examples/extended/persistency/gdml/G01	      	example-ext-persistency-*	   g01.out
#update  examples/extended/persistency/gdml/G02	      	example-ext-persistency-*	   read_gdml.out
#update  examples/extended/persistency/gdml/G03	      	example-ext-persistency-*	   read_ext.out
update  examples/extended/persistency/gdml/G04	      	example-ext-persistency-*	   g04.out
update  examples/extended/persistency/P03	      	example-ext-persistency-*	   batch.out

update  examples/extended/radioactivedecay/rdecay01	example-ext-radioactivedecay-*      rdecay01.out
update  examples/extended/radioactivedecay/rdecay02	example-ext-radioactivedecay-*      rdecay02.out

update  examples/extended/runAndEvent/RE01		example-ext-runandevent-*      sample.out
update  examples/extended/runAndEvent/RE02		example-ext-runandevent-*	   run.out
#update  examples/extended/runAndEvent/RE02		*run3-run  run3.out
#update  examples/extended/runAndEvent/RE02		*run4-run  run4.out
update  examples/extended/runAndEvent/RE03		example-ext-runandevent-*	   run.out

#update  examples/extended/visualization/userVisAction	*run1-run	run1.out
#update  examples/extended/visualization/userVisAction	*run2-run	run2.out
#update  examples/extended/visualization/userVisAction	*uva-run	userVisAction.out

update  examples/advanced/amsEcal	     		example-adv-*	   run1.out
update  examples/advanced/dnaphysics	     		example-adv-*	   dna.out
#update  examples/advanced/eRosita	      		example-adv-*	   eRosita.out
update  examples/advanced/gammaray_telescope   		example-adv-*	   gammaraytel.out
update  examples/advanced/microbeam      		example-adv-*	   microbeam.out
update  examples/advanced/microdosimetry		example-adv-*	   microdosimetry.out
update  examples/advanced/nanobeam      		example-adv-*	   nanobeam.out
#update  examples/advanced/cosmicray_charging      	example-adv-*	   shoot.out

#basic examples
update  examples/basic/B1      				example-bas-*	   exampleB1.out
update  examples/basic/B2/B2a  				example-bas-*	   exampleB2a.out
update  examples/basic/B2/B2b  				example-bas-*	   exampleB2b.out
update  examples/basic/B3      				example-bas-*	   exampleB3.out
update  examples/basic/B4/B4a  				example-bas-*	   exampleB4a.out
update  examples/basic/B4/B4b  				example-bas-*	   exampleB4b.out
update  examples/basic/B4/B4c  				example-bas-*	   exampleB4c.out
update  examples/basic/B4/B4d  				example-bas-*	   exampleB4d.out
