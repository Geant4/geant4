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
print "/vis/open/VRML1" ;
print "/vis/draw/current" ;
print "/vis~/create_view/new_graphics_system VRML1" ;
print "/test/errorFileName " ARGV[1] ;
}

# print the solid command
NR == 2 { print substr($0, 3) ;}
{ print "#" $0 ; } 
$1 ~ /[0-9]+/ {
 print "/test/draw " $1 ;
 print "/test/pause" ;
}
