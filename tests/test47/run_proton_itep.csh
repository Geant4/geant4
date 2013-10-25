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

set lists = ( "bertini" "binary" "ftfp" )

set energyPRO = ( 1400 7500 )
set targetPRO = ( "C" "U" )
set energyPRO_GeV = ( 1.40 7.50 )


set NEnTotal = $#energyPRO
set NTgTotal = $#targetPRO

if ( $#argv > 0 ) then

set ii=$1
@ ii = ( $ii + 1 )

if ( $ii > 4 ) then
@ ii = ( $ii - 8 )
endif

@ i = ( $ii + 1 )
@ i = ( $i / $NTgTotal )
@ NE = $i
@ jj = ( $ii + 1 )
@ jj = ( $jj % $NTgTotal )
@ jj = ( $jj + 1 )
@ NT = $jj

else

set i  = 1
set NE = $#energyPRO
set jj = 1
set NT = $#targetPRO

endif

while ( $i <= $NE )

set j=$jj
while ( $j <= $NT )

cat > ITEP.pro.$energyPRO[$i].$targetPRO[$j] <<EOF
#verbose
0
#rad
#events
1000000
//--------Proton_processes
#particle
proton
//--------
#isITEP
#position(mm)
0. 0. 0.
#direction
0. 0. 1.
//--------
#material
$targetPRO[$j]
#momentum(MeV/c)
$energyPRO[$i]
// ---
EOF

if ( $#argv > 1 )
{
set seed=1234
@ seed = ( $seed + $2 )
printf "myseed\n" >> ITEP.pro.$energyPRO[$i].$targetPRO[$j]
printf "%d\n" $seed >> ITEP.pro.$energyPRO[$i].$targetPRO[$j]
printf "#jobID\n" >> ITEP.pro.$energyPRO[$i].$targetPRO[$j]
printf "%d\n" $2 >> ITEP.pro.$energyPRO[$i].$targetPRO[$j]
}

if ( $energyPRO[$i] == 1400 ) then
set maxLists = 3
else
set maxLists = 5
endif

set k=0
while ( $k != $maxLists )
@ k = $k + 1
printf "#generator\n" >> ITEP.pro.$energyPRO[$i].$targetPRO[$j]
printf "%s\n" $lists[$k] >> ITEP.pro.$energyPRO[$i].$targetPRO[$j]
printf "#run\n" >> ITEP.pro.$energyPRO[$i].$targetPRO[$j]
end 
printf "#exit\n" >> ITEP.pro.$energyPRO[$i].$targetPRO[$j]

$G4EXE/test47 ITEP.pro.$energyPRO[$i].$targetPRO[$j]

@ index = ( $i - 1 ) 
@ index = ( $index * $NTgTotal )
@ index = ( $index  + $j )
@ index = ( $index + 8 )

$ROOTSYS/bin/root -b -p -q Plot.C\($index\) 

set outputPlots = proton$targetPRO[$j]toprotonat$energyPRO_GeV[$i]GeV_1.eps 

### these actions below are just examples of moving the output around...
###
###/usr/bin/X11/gv $outputPlots &
###/bin/cp $outputPlots /uscms/home/yarba_j/G4test47/.


if ( -e ITEP.pro.$energyPRO[$i].$targetPRO[$j] ) then
rm ITEP.pro.$energyPRO[$i].$targetPRO[$j]
endif


@ j = $j + 1
end
@ i = $i + 1
end
