#!/bin/tcsh -f

# this is PROTOTYPE script
# at present, it assumes that all G4 (and related) enviroment variables 
# are set manually, for example via setup [tool] [version] [...] command,
# plus G4EXE env.variable is also set;
# however, in order to execute in batch/grid, setup of G4 (etc.) enviroment
# must be done EXPLICITLY

# the following settings have been done to run in the OS Grid,
# submitted via CMSLPC cluster;
# application has been composed against CMS build of G4.9.2
#
###setenv G4INSTALL /uscmst1/prod/sw/cms/slc4_ia32_gcc345/external/geant4/9.2
###setenv ROOTSYS /uscmst1/prod/sw/cms/slc4_ia32_gcc345/lcg/root/5.22.00a-cms4/
###setenv CLHEP_BASE_DIR /uscmst1/prod/sw/cms/slc4_ia32_gcc345/external/clhep/1.9.4.2

#
# this sets it up to CMS build of G4 on LXPLUS
#
setenv G4INSTALL /afs/cern.ch/cms/sw/slc4_ia32_gcc345/external/geant4/9.2
setenv ROOTSYS /afs/cern.ch/cms/sw/slc4_ia32_gcc345/lcg/root/5.22.00a-cms5
setenv CLHEP_BASE_DIR /afs/cern.ch/cms/sw/slc4_ia32_gcc345/external/clhep/1.9.4.2

setenv CLHEP_LIB CLHEP

setenv G4SYSTEM Linux-g++
setenv G4LIB    $G4INSTALL/lib
setenv G4DATA   $G4INSTALL/data

setenv G4LEVELGAMMADATA $G4DATA/PhotonEvaporation2.0
setenv G4NEUTRONHPDATA $G4DATA/G4NDL3.12
setenv G4RADIOACTIVEDATA $G4DATA/RadioactiveDecay3.2
setenv G4LEDATA $G4DATA/G4EMLOW5.1

# this is for CMSLPC
#
###cd /uscms_data/d2/yarba_j/geant4/tests
#
# this is for LXPLUS
#
cd $HOME/scratch0/geant4/tests

setenv G4WORKDIR  $PWD
setenv G4EXE $G4WORKDIR/bin/$G4SYSTEM

setenv LD_LIBRARY_PATH $G4LIB/:$CLHEP_BASE_DIR/lib/:$ROOTSYS/lib/

cd test47 

set lists = ( "bertini" "binary" "ftfp"  )

set energyPIM = 5000 
set targetPIM = ( "C" "Cu" "Pb" "U" )
set energyPIM_GeV = 5.00 


set NTgTotal = $#targetPIM

if ( $#argv > 0 ) then

set ii=$1
@ ii = ( $ii + 1 )

set i = $ii

if ( $i > 4 ) then
@ i = ( $i - 4 )
endif

set NT = $i

else

set i = 1
set NT = $#targetPIM

endif

set maxLists = $#lists

while ( $i <= $NT )

cat > ITEP.pim.$energyPIM.$targetPIM[$i] <<EOF
#verbose
0
#rad
#events
1000000
//--------Proton_processes
#particle
pi-
//--------
#isITEP
#position(mm)
0. 0. 0.
#direction
0. 0. 1.
//--------
#material
$targetPIM[$i]
#momentum(MeV/c)
$energyPIM
// ---
EOF

set k=0
while ( $k != $maxLists )
@ k = $k + 1
printf "#generator\n" >> ITEP.pim.$energyPIM.$targetPIM[$i]
printf "%s\n" $lists[$k] >> ITEP.pim.$energyPIM.$targetPIM[$i]
printf "#run\n" >> ITEP.pim.$energyPIM.$targetPIM[$i]
end 
printf "#exit\n" >> ITEP.pim.$energyPIM.$targetPIM[$i]

$G4EXE/test47 ITEP.pim.$energyPIM.$targetPIM[$i]

@ index = ( $i + 4 )

$ROOTSYS/bin/root -b -p -q Plot.C\($index\) 

set extension = GeV_1.eps
set outputPlots1 = piminus$targetPIM[$i]toprotonat$energyPIM_GeV$extension
set outputPlots2 = piminus$targetPIM[$i]toneutronat$energyPIM_GeV$extension

### these actions below are just examples of moving the output around...
###
###/usr/bin/X11/gv $outputPlots1 &
###if ( -e $outputPlots1 ) /bin/cp $outputPlots1 /uscms/home/yarba_j/G4test47/.
###if ( -e $outputPlots2 ) /bin/cp $outputPlots2 /uscms/home/yarba_j/G4test47/.


if ( -e ITEP.pim.$energyPIM.$targetPIM[$i] ) then
rm ITEP.pim.$energyPIM.$targetPIM[$i]
endif

@ i = $i + 1
end
