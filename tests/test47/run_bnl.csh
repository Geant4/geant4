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

set lists = ( "bertini" "qgsp" "ftfp" )

set energy = 14600 
set target = ( "Be" "Cu" "Au" )
set energy_GeV = 14.6 


set NTgTotal = $#target

if ( $#argv > 0 ) then

set ii=$1
@ ii = ( $ii + 1 )

set i = $ii

if ( $i > 3 ) then
@ i = ( $i - 12 )
endif

set NT = $i

else

set i = 1
set NT = $#target

endif

set maxLists = $#lists

while ( $i <= $NT )

cat > BNL.proton.$energy.$target[$i] <<EOF
#verbose
0
#rad
#events
1000000
//--------Proton_processes
#particle
proton
//--------
#isBNL
#position(mm)
0. 0. 0.
#direction
0. 0. 1.
//--------
#material
$target[$i]
#momentum(MeV/c)
$energy
// ---
EOF

set k=0
while ( $k != $maxLists )
@ k = $k + 1
printf "#generator\n" >> BNL.proton.$energy.$target[$i]
printf "%s\n" $lists[$k] >> BNL.proton.$energy.$target[$i]
printf "#run\n" >> BNL.proton.$energy.$target[$i]
end 
printf "#exit\n" >> BNL.proton.$energy.$target[$i]

$G4EXE/test47 BNL.proton.$energy.$target[$i]


$ROOTSYS/bin/root -b -p -q Plot.C\($i,\"BNL\"\) 

set extension = GeV.eps
set outputPlots1 = proton$target[$i]topiplusat$energy_GeV$extension
set outputPlots2 = proton$target[$i]topiminusat$energy_GeV$extension
set outputPlots3 = proton$target[$i]toprotonat$energy_GeV$extension
set outputPlots4 = proton$target[$i]tokplusat$energy_GeV$extension
set outputPlots5 = proton$target[$i]tokminusat$energy_GeV$extension

### these actions below are just examples of moving the output around...
###
###/usr/bin/X11/gv $outputPlots1 &
###if ( -e outputPlots1) /bin/cp $outputPlots1 /uscms/home/yarba_j/G4test47/.


if ( -e BNL.proton.$energy.$target[$i] ) then
rm BNL.proton.$energy.$target[$i]
endif

@ i = $i + 1
end
