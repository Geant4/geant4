#!/bin/tcsh -f

source g4setup_grid.csh

set lists = ( "lepar" "bertini" "binary" "ftfp" "qgsc" )

# defaults...
set Beam = ( "proton" )
set Energy = ( 1400 )
set Target = ( "C" )
set NEvents = ( 10000 )
set JobID = ( -1 )
set ClusterID = ( -1 )


# settings
if ( $#argv > 0 ) then 
set Beam = $1
set Energy = $2
set Target = $3
set NEvents = $4
set JobID = $5
set ClusterID = $6
endif
cat > ITEP.$Beam.$Energy.$Target.$JobID.$ClusterID <<EOF
#verbose
0
#rad
#events
$NEvents
//--------Proton_processes
#particle
$Beam
//--------
#isITEP
#position(mm)
0. 0. 0.
#direction
0. 0. 1.
//--------
#material
$Target
#momentum(MeV/c)
$Energy
// ---
EOF

set seed=1234
@ seed = ( $seed + $JobID )
printf "myseed\n" >> ITEP.$Beam.$Energy.$Target.$JobID.$ClusterID
printf "%d\n" $seed >> ITEP.$Beam.$Energy.$Target.$JobID.$ClusterID
printf "#jobID\n" >> ITEP.$Beam.$Energy.$Target.$JobID.$ClusterID
printf "%d\n" $JobID >> ITEP.$Beam.$Energy.$Target.$JobID.$ClusterID
printf "#clusterID\n" >> ITEP.$Beam.$Energy.$Target.$JobID.$ClusterID
printf "%d\n" $ClusterID >> ITEP.$Beam.$Energy.$Target.$JobID.$ClusterID

if ( $Energy == 1400 ) then
set maxLists = 3
else
set maxLists = 5
endif

set k=0
while ( $k != $maxLists )
@ k = $k + 1
printf "#generator\n" >> ITEP.$Beam.$Energy.$Target.$JobID.$ClusterID
printf "%s\n" $lists[$k] >> ITEP.$Beam.$Energy.$Target.$JobID.$ClusterID
printf "#run\n" >> ITEP.$Beam.$Energy.$Target.$JobID.$ClusterID
end 
printf "#exit\n" >> ITEP.$Beam.$Energy.$Target.$JobID.$ClusterID

$G4EXE/test47 ITEP.$Beam.$Energy.$Target.$JobID.$ClusterID

if ( -e ITEP.$Beam.$Energy.$Target.$JobID.$ClusterID ) then
###rm ITEP.$Beam.$Energy.$Target.$ID
endif


