#!/usr/bin/perl
#
#  Process the output of a run in which checking for Energy and momentum conservation is activated.
#  Summarise each issue reported in one of two ways : 
#      either from default check for 'catastrophic' violation of Energy conservation
#      or for problem with 'developer' level checks of E/p, charge and baryon number conservation.
#
#  NOTE: Since 9.6-beta both relative and absolute E/p checks must fail for an event to report a failure
#        (  Version before 9.5-ref-04 or so reported an issues if only one failed. )
#
#  First   version:   John Apostolakis,  31 May  2012   ( in AWK ) 
#  Latest revision:   John Apostolakis,  29 June 2012
#

#  INPUT Format:
#rint process, model, primary, primaryEn, nucleus, limitR, energyR, momentumR, limitA, energyA, momentumA, chargeB, baryonB; 
#
#   1        2              3      4      5            6      7       8        9   10       11  12        13    14
#  Process: ElectroNuclear  ,      Model: CHIPSElectroNuclear
#  Primary: e-            (11),    E=     39.5593,     target nucleus (82,207)
#  fail     relative,     limit    0.01,  values       E      /       p     (MeV)  =  0.390098  /   0.385144
#  fail     absolute,     limit    (MeV)  1,           values E       /         p  MeV)      =  5.2327     /   15.2347
#  pass     charge/baryon number   balance 0 / 0 

eval 'exec /usr/bin/perl -S $0 ${1+"$@"}'
    if $running_under_some_shell;
			# this emulates #! processing on NIH machines.
			# (remove #! line above if indigestible)

eval '$'.$1.'$2;' while $ARGV[0] =~ /^([A-Za-z_0-9]+=)(.*)/ && shift;
			# process any FOO=bar switches


$[ = 1;			# set array base to 1
$, = ' ';		# set output field separator
$\ = "\n";		# set output record separator

$verbose = 1;

printf "%-18s %-20s %6s  %8s  %9s, %8s  %9s %9s %8s  %9s %9s  %3s  %3s\n", 
'Process', 'Model', 'primary', 'primaryEn', 'nucleus', 'limitR', 'energyR',

  'momentumR', 'limitA', 'energyA', 'momentumA', 'chargeB', 'baryonB';

## Original lines:
#  Process: ElectroNuclear , Model: CHIPSElectroNuclear
# Primary: e- (11), E= 39.5593, target nucleus (82,207)
#   fail relative, limit 0.01, values E / p (MeV) = 0.390098 / 0.385144
#   fail absolute, limit (MeV) 1, values E / p (MeV) = 15.2327 / 15.2347
#   pass charge/baryon number balance 0 / 0 

while (<>) {
    chomp;	# strip record separator
    @Fld = split(' ', $_, -1);
    if (/Process/ && /Model/) {
	$process = $Fld[2];
	$model = $Fld[5];
    }
    if (/Primary/) {
	$primary = $Fld[2];
	$primaryEn = $Fld[5] =~ s/,//;
	$nucleus = $Fld[8];
    }
    ##            (11), E= 39.5593, target nucleus (82,207)
    if (/relative/) {
	$limitR = $Fld[4];
	$energyR = $Fld[11];
	$momentumR = $Fld[13];
    }
    if (/absolute/) {
	$limitA = $Fld[5];
	$energyA = $Fld[12];
	$momentumA = $Fld[14];
	$energyR = $deltaE / $primaryEn;
    }
    if (/charge/) {
	$chargeB = $Fld[5];
	$baryonB = $Fld[7];
    }
    if (/charge/) {
	&report_it();
    }

    ## Problem is found in Excpetion
    if (/G4Exception/ && /START/) {
	$InException = 1;
    }

    ## /Primary/ && InException { particle= $2; pdg= $11; En=$5;  nucleus=$8; } 
    if (/initial/ && /final/ && $InException) {
	$deltaE = $Fld[5];
	$unit = $Fld[6];

	$energyA = $deltaE;

	$energyR = $deltaE / $primaryEn;

	$momentumA = $momentumR = 0.0;

	$chargeB = 0;
	$baryonB = 0;
    }

    ##  E(initial - final) = 51787.2 MeV.
    if (/Process/ && /Model/ && $InException) {
	$process = $Fld[4];
	$model = $Fld[6];

	$limitR = 0;
	$limitA = 5000;
    }
    ##  Process / Model: ElectroNuclear / CHIPSElectroNuclear

    if (/G4Exception/ && /END/) {
	## print "*** From Exception *** "; 
	&report_it();

	$InException = 0;
    }

    ## E(initial - final) = 15001.8 MeV.

    ## -------- WWWW ------- G4Exception-START -------- WWWW -------
    ## *** G4Exception : had012
    ##       issued by : G4HadronicProcess:CheckResult()
    ## Warning: Bad energy non-conservation detected, will re-sample the interaction
    ##  Process / Model: ElectroNuclear / CHIPSElectroNuclear
    ##  Primary: e- (11), E= 44754.6, target nucleus (82,208)
    ##  E(initial - final) = 15001.8 MeV.
    ## 
    ## *** This is just a warning message. ***
    ## -------- WWWW -------- G4Exception-END --------- WWWW -------

    if (/Event/) {
	if ($verbose) {
	    print $_;
	}
    }
    if (/\/gun/ || /\gun\/energy/) {
	if ($verbose) {
	    print $_;
	}
    }
    if (/Material/) {
	if ($verbose) {
	    print $_;
	}
    }
}

sub report_it() {
	#print process, model, primary, primaryEn, nucleus, limitR, energyR, momentumR,  limitA, energyA, momentumA, chargeB, baryonB; 
	#print prox  mod  pr  prEn  nucl  limR enR,  pR,  limA  enA   momA  chargeB baryonB
	#printf   "%-18s %-20s %6s  %8.2f  %9s, %8s  %9.4f ", ## " %9.4f %8s  %9.2f %9.2f  %3s  %3s\n",
	#          $process, $model, $primary, $primaryEn, $nucleus, $limitR, $energyR;
	printf
	  "%-18s %-20s %6s  %8.2f  %9s, %8s  %9.4f  %9.4f %8s  %9.4f %9.4f  %3s  %3s\n",
	   $process, $model, $primary, $primaryEn, $nucleus, $limitR, $energyR,
	   $momentumR, $limitA, $energyA, $momentumA, $chargeB, $baryonB;
        # print $limitR;
        # print $energyR;
        # print $momentumR;
        # print $limitA;
        # printf  "%9.4f %8s  %9.2f %9.2f  %3s  %3s\n",
	#    $momentumR, $limitA, $energyA, $momentumA, $chargeB, $baryonB;

	#rint prox mod           pr     prEn     nucl     limR  energyR, pR       lA energyA limA, chargeB, baryonB; 
}
# ElectroNuclear CHIPSElectroNuclear  30.5987, (82,208) 0.01, 0.352011 0.346255 1, 10.5912 10.5935 0 0

