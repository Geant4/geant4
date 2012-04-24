#!/usr/bin/perl -w 

# Author: Oliver Link
#
# Mon May 23 18:02:09 CEST 2005
#
# Script to test systematically solids with SurfaceChecker
# For every solid a hbook file will be generated in folder hbk
#   eg. for 500 events , G4Box:   hbk/run_Box_500.hbk
# Attention: already existing hbk files will be replaced

use strict ;

# list of solids for testing
my @solids = qw/Tet Trap Torus Box Sphere Tubs Orb Cons TwistedTubs TwistedBox
    TwistedTrd TwistedTrap  Ellipsoid EllipticalCone EllipticalTube Hype Shell
    HalfSphere HollowSphere b1Ub2 b1Ib2 b1Sb2 b1Ub1 b1Ib1 b1Sb1/ ;

@solids = @ARGV if (($#ARGV+1)>0);	# Process solids specified by user

my $nevents = 100000 ;   # sets the number of events

# --------------------------------------------------------------------------


my $macro ;    # the name of the macro file.    Will be deleted at the end.
my $data ;     # output file of SurfaceChecker. Will be deleted at the end.
my $hbk  ;     # resulting hbook file.

# Ensure that data/ and hbk/ directories exist
unless (-d "data") { mkdir "data" or die; }
unless (-d "hbk")  { mkdir "hbk" or die; }

foreach my $solid ( @solids ) {

    $macro = "run_$solid" . ".mac" ;
    $data  = "data/run_$solid" . "_$nevents" . ".data" ;
    $hbk   = "hbk/run_$solid"  . "_$nevents" . ".hbk" ;

    $data =~ tr/A-Z/a-z/ ;     # paw dosn't like capital letters for files
    $hbk  =~ tr/A-Z/a-z/ ;     
 
    print "process solid $solid: $macro  --> $data\n" ;

#   prepare the macro file for the G4 application

    open(FH,">$macro") || die "cannot open file $macro:$!\n" ;

    print FH "/run/initialize\n" ;
    print FH "/mydet/SelectDetector $solid\n" ;
    print FH "/gun/particle geantino\n" ;
    print FH "/run/beamOn $nevents\n" ;

    close(FH) || die "cannot close file $macro:$!\n" ;

#   run SurfaceChecker application and redirect output to data file
    
    system("SurfaceChecker $macro > $data") && die "cannot execute SurfaceChecker for macro $macro to file $data:$!\n" ;

#   Convert data to hbook file. Attention: already existing hbk files will be deleted
    if ( -e $hbk ) {
	print " >>>> Attention: file $hbk exists already. It will be replaced by a new version...\n" ;
	unlink($hbk) || die "cannot delete file $hbk:$!\n" ;
    }
    system("extract.pl > /dev/null") && die "cannot execute extract.pl:$!\n" ;

# delete temporary macro/data files. Uncomment the lines if you want to keep the files.

    unlink($macro) || die "cannot delete file $macro:$!\n" ;
    unlink($data)  || die "cannot delete file $data:$!\n" ;

}

system(" paw -w 0 -b solid > /dev/null") && die "cannot execute paw macro solid:$!\n" ;
