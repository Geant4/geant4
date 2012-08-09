#!/usr/bin/perl
eval 'exec /usr/bin/perl -S $0 ${1+"$@"}'
    if $running_under_some_shell;
			# this emulates #! processing on NIH machines.
			# (remove #! line above if indigestible)

eval '$'.$1.'$2;' while $ARGV[0] =~ /^([A-Za-z_0-9]+=)(.*)/ && shift;
			# process any FOO=bar switches

#
#  Given output of test22 with '/tracking/verbose 1', produce a table of the
#  number of times each interaction occurs.
#  Displays the Material, primary particle and Energy (value, unit), 
#   in addition to output data.
#
#  Run it as follows: 
#     $G4BIN/Linux-g++/test22 test22.it10-verbose.in | ./process-GammaN.pl
#
# First version: John Apostolakis, 7 June 2012
#
$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator

$started = 0;
$verbose = 1;
## Not started

$debug = 1;    ## 

$Primary="None"; 
$Energy=0.0; 
$Unit="None";
$Events= $numPhotonA= $numElectronA= $numPositronA= 0;
$errPhotonA= $errElectronA= $errPositronA = 0;

## &banner();
printf "%10s  %10s  %7s  %5s  %8s :  %8s %8s %8s    %8s %8s %8s \n", 'Material',
  'Primary', 'Energy', 'Unit', 'Events', 'PhotonA', 'ElectronA', 'PositronA', 'errPhotonA', 'errElectronA', 'errPositronA';


while (<>) {
    ($Fld1,$Fld2,$Fld3,$Fld4,$Fld5,$Fld6,$Fld7,$Fld8,$Fld9,$Fld10) = split(' ', $_, -1);
    if (/PhotonInelastic/) {
        if( /Process/ ){ 
	   $errPhotonA++;
        }else{ 
	   $numPhotonA++;
        }

	if ($verbose > 6) {
	    print $Fld10;
	}
    }
    if (/ElectroNuclear/) {
        if( /Process/ ){ 
	   $errElectronA++;
        }else{ 
	   $numElectronA++;
	}
    }
    if (/PositronNuclear/) {
        if( /Process/ ){ 
	    $errPositronA++;
	}else{
	    $numPositronA++;
	}
    }

    if (m#/gun# || m#/run# || m#/SelectMaterial# ) {
	if ($started) {
	    &report_last_energy();
	    ## &clear_sums();
	    $numPhotonA = 0;
	    $numElectronA = 0;
	    $numPositronA = 0;
	    $errPhotonA = 0;
	    $errElectronA = 0;
	    $errPositronA = 0;

	    $started = 0;
	}
    }

    if (/beamOn/ ) {
	if (/Idle/) {
	    $Events = $Fld3;
	}
	else {
	    $Events = $Fld2;
	}
	$started = 1;

	if ($debug) {
	    printf "Dbg> Line: %9d   Beam On   Events=%8d\n", $., $Events;
	}
    }

    if (/gun/ && /energy/) {
	if (/Idle/) {
	    $Energy = $Fld3;
	    $Unit = $Fld4;
	}
	else {
	    $Energy = $Fld2;
	    $Unit = $Fld3;
	}
	if ($debug) {
	    printf "Dbg> Line: %9d   Energy=         %10s   Unit=%8s\n", $., $Energy,

	      $Unit;
	}
    }

    if (/SelectMaterial/) {
        if ( ! /Idle/ ) { 
	    $Material = $Fld2;
	}else{
	    $Material = $Fld3;
	}
	if ($debug) {
	    printf "Dbg> Line: %9d   Select Material= %10s \n", $., $Material;
	}
    }

    if (/gun\/particle/) {
        if ( ! /Idle/ ) { 
	    $Primary = $Fld2;
	}else{
	    $Primary = $Fld3;
	}

	if ($debug) {
	    printf "Dbg> Line: %9d   Particle=        %10s \n", $., $Primary;
	}
        # print "Material = ", $Primary    if $debug;
        # print "Line>  ", $_              if $debug;
    }

    # /run/verbose 2
    # /run/initialize
    # Bias factor for Phot Nuclear is set to 100
    # Bias factor for ElectronNuclear is set to 1000
    # Bias factor for PositronNuclear is set to 1000 - same process instance is used as for electron 
    # /gun/direction 1 0 0
    # /tracking/verbose 1
    # /mydet/SelectMaterial Pb
    # /gun/particle gamma
    # 
    # /gun/energy 1.5 GeV
    # /run/beamOn 10
    #      ElectroNuclear  Models:          CHIPSElectroNuclear: Emin(GeV)=    0  Emax(GeV)= 30000
    #      ElectroNuclear  Crs sctns:          ElectroNuclearXS: Emin(GeV)=    0  Emax(GeV)= 100000
    #      PhotonInelastic  Models:            CHIPSGammaNuclear: Emin(GeV)=    0  Emax(GeV)= 3.5
    #      PhotonInelastic  Crs sctns:            PhotoNuclearXS: Emin(GeV)=    0  Emax(GeV)= 100000
    # Event number 0
    #      1      199        0        0         0        0      199       199       SPhys PhotonInelastic
}
&report_last_energy();

sub report_last_energy() {
    if( $debug ) { printf  "Dbg> %8d : Called Report_Last_Energy() \n", $.;  } 
    printf "%10s  %10s  %7g  %5s  %8d :  %8d %8d %8d   %8d %8d %8d \n", $Material,
      $Primary, $Energy, $Unit, $Events, $numPhotonA, $numElectronA, $numPositronA,  
					 $errPhotonA, $errElectronA, $errPositronA  ;

    printf "XLout> %10s,  %10s,  %7g,  %5s,  %8d, :  %8d, %8d, %8d,   %8d, %8d, %8d \n", $Material,
      $Primary, $Energy, $Unit, $Events, $numPhotonA, $numElectronA, $numPositronA,  
					 $errPhotonA, $errElectronA, $errPositronA  ;
}
