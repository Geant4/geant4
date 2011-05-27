#!/bin/csh -f
# Main file for production 

@ E = 0
set Cut = ( 1 0 )
#set L = ( 1 10 20 50 100 )
set L = ( 20 80 )
#set F = ( QGSP_BERT FTFP_BERT LHEP CHIPS QBBC QGSP_BIC QGSP_INCL_ABLA Shielding )
set F = ( QGSP_BERT FTFP_BERT LHEP CHIPS QGSP_BIC )
set Z = ( 2 4 6 8 10 14 16 18 26 54 82 92 ) #12
set A = ( 4 9 12 16 20 28 32 40 56 131 207 238 )
#set e = ( 10 50 100 150 200 250 350 500 800 1000 2000 3000 4000 7000 10000 ) #14
set e = ( 10 60 250 500 800 1200 3800 10000 ) #8
#set Mat = ( He B C O Al Ar Fe Cu Sm Ta Pb U WATER CESIUM_IODIDE POLYETHYLENE MYLAR SILICON_DIOXIDE KEVLAR STAINLESS-STEEL AIR ) #20
set Mat = ( He B C O Al Ar Ti Fe Cu Zr Sm Ta Pb U WATER CESIUM_IODIDE POLYETHYLENE SILICON_DIOXIDE KEVLAR AIR ) #20


set T = 2 
while ($T < 3)
set l = 1 
while ($l < 3)
setenv REFERENCE geant4-09-04-ref-03m_$Cut[$T]_$L[$l]
mkdir -p $RESULTS/HADR02/$REFERENCE/ 
set i = 2
while ($i < 3)
setenv PHYSLIST $F[$i]
echo "PRODUCTION in " $RESULTS/HADR02/$REFERENCE/$F[$i]
mkdir -p $RESULTS/HADR02/$REFERENCE/$F[$i]/
echo $F[$i]

set n = 1 
while ($n < 9)
set m = 1 
while ($m < 21)
set j = 2 
while ($j < 13)
@ E = ( $e[$n] * $A[$j] )

cd $HADR02
echo "/control/verbose 1" >! run.in
echo "/run/verbose 1" >> run.in
echo "/tracking/verbose 0" >> run.in
echo "/testhadr/TargetMat " G4_$Mat[$m]  >> run.in
echo "/testhadr/TargetRadius  20 cm"  >> run.in
echo "/testhadr/TargetLength  " $L[$l] " cm"  >> run.in
echo "/testhadr/PrintModulo   1000"  >> run.in
echo "/testhadr/ionPhysics DPMJET"  >> run.in
echo "/testhadr/HistoName"  $A[$j]_$Mat[$m]_$E  >> run.in
echo "/run/initialize"  >> run.in
if ($Cut[$T] > 0 ) then 
echo "/run/setCut"   $Cut[$T] " km"  >> run.in
endif
echo "/gun/particle ion"  >> run.in
echo "/gun/ion " $Z[$j] " " $A[$j]  >> run.in
echo "/gun/energy " $E " GeV"  >> run.in
echo "/run/beamOn 2000 "  >> run.in

echo Material $Mat[$m] , Energy $E GeV , beam ion $Z[$j] $A[$j] , target length $L[$l] cm , $F[$i]
cat run.in
cd $RESULTS/HADR02/$REFERENCE/$F[$i]/
$G4WORKDIR/bin/Linux-g++/hadr02 $HADR02/run.in >&! $Mat[$m]_$Z[$j]_$E.out 
echo
echo
echo
echo

@ j++
end
@ m++
end
@ n++
end
@ i++
end 
@ l++
end 
@ T++
end 
exit





