#!/bin/awk -f
# MEDERNACH Emmanuel Aug 2000
#
# This script create a script for SBT to 
# show errors. You need to provide the solid
#
#
# Usage : ErrorView <log file> > file.SBT

BEGIN { 
print "# ErrorView generated script for viewing error in SBT #" ;
print "/vis/open VRML2FILE"
print "/vis/viewer/set/style wireframe"
print "#/vis/open VRML1FILE"
print "/vis/drawVolume test"
print "/tracking/storeTrajectory 1"
print "/test/errorFileName " ARGV[1] ;
}

# print the solid command
NR == 2 { print substr($0, 3) ;}
{ print "#" $0 ; } 
$1 ~ /[0-9]+/ {
 print "/test/draw " $1 ;
 print "/test/pause" ;
}
