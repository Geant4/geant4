#!/usr/local/bin/perl

# stt.geant4-05-00-cand-00+tags1023
#$SttDir2 = "$ENV{'G4WORKDIR'}/stt.geant4-05-00-cand-00+tags1023/$ENV{'G4SYSTEM'}/";
# ref-04
#$SttDir2 = "$ENV{'G4WORKDIR'}/stt.geant4-04-01-ref-03+tags995/$ENV{'G4SYSTEM'}/";
# 
#$SttDir2 = "$ENV{'G4WORKDIR'}/stt.geant4-05-00-cand-01+tags1049/$ENV{'G4SYSTEM'}/";
#$SttDir2 = "$ENV{'G4WORKDIR'}/stt.geant4-05-00-cand-00+tags1021/$ENV{'G4SYSTEM'}/";
#$SttDir2 = "$ENV{'G4WORKDIR'}/stt.geant4-05-00-cand-00+tags1024/$ENV{'G4SYSTEM'}/";
#$SttDir = "$ENV{'G4WORKDIR'}/stt/$ENV{'G4SYSTEM'}/";
#$SttDir = "$ENV{'G4WORKDIR'}/stt.geant4-04-01-ref-04+tags1019//$ENV{'G4SYSTEM'}/";
#$SttDir2 = "$ENV{'G4WORKDIR'}/stt.geant4-04-01-ref-04+tags1019//$ENV{'G4SYSTEM'}/";
#    $SttDir2 = "$ENV{'G4WORKDIR'}/stt.geant4-05-00-cand-01+tags1049/$ENV{'G4SYSTEM'}/$TestName";

#    $SttDir2 = "/afs/cern.ch/sw/geant4/stt/dev1/Linux-g++/optim_7.3_3.2/stt.geant4-05-00-ref-03+tags1728/$ENV{'G4SYSTEM'}";

    $SttDir2 = "/afs/cern.ch/sw/geant4/stt/dev1/Linux-g++/optim_7.3_3.2/stt.geant4-05-00-ref-03+tags1799/$ENV{'G4SYSTEM'}";


#    $SttDir2 = "$ENV{'G4WORKDIR'}/stt.geant4-05-00-cand-01+tags1062/$ENV{'G4SYSTEM'}/$TestName";


#$SttDir2 = $SttDir;
#$SttDir2 =~ s|prod|dev1|;

#print "SttDir  = $SttDir\n";
#print "SttDir2 = $SttDir2\n";

my $deltatime1 ;
my $deltatime2 ;

opendir(TT,"$SttDir2") || die "Failed to opendir  $SttDir2 $!";
@Tests = grep(m/\.err$/,readdir(TT));
closedir(TT);

foreach $Test (@Tests) {
#

$TestName = $Test;

print "TestName = $TestName\n";

#
# Get the size of $TestName.err 
#
#stt.geant4-05-00-cand-00+tags1021
#stt.geant4-05-00-cand-01+tags1049/
#stt.geant4-04-01-ref-05+tags1073/

#    $filename = "$ENV{'G4WORKDIR'}/stt.geant4-05-00-cand-01+tags1049//$ENV{'G4SYSTEM'}/$TestName";
#    $filename = "$ENV{'G4WORKDIR'}/stt.geant4-05-00-cand-00+tags1031/$ENV{'G4SYSTEM'}/$TestName";
#    $filename = "$ENV{'G4WORKDIR'}/stt.geant4-04-01-ref-04+tags1005/$ENV{'G4SYSTEM'}/$TestName";
#    $filename = "$ENV{'G4WORKDIR'}/stt.geant4-05-00-cand-00+tags1021/$ENV{'G4SYSTEM'}/$TestName";
#    $filename = "$ENV{'G4WORKDIR'}/stt.geant4-05-00-cand-01+tags1041/$ENV{'G4SYSTEM'}/$TestName";
#    $filename = "$ENV{'G4WORKDIR'}/stt.geant4-05-00-cand-01+tags1053//$ENV{'G4SYSTEM'}/$TestName";
#    $filename = "$ENV{'G4WORKDIR'}/stt.geant4-05-00-cand-01+tags1062/$ENV{'G4SYSTEM'}/$TestName";
#    $filename = "$ENV{'G4WORKDIR'}/stt.geant4-04-01-ref-05+tags1073/$ENV{'G4SYSTEM'}/$TestName";

#    $filename = "/afs/cern.ch/sw/geant4/stt/dev2/Linux-g++/optim_7.3_3.2/stt.geant4-05-00-ref-03+tags1729/$ENV{'G4SYSTEM'}/$TestName";

    $filename = "/afs/cern.ch/sw/geant4/stt/dev1/Linux-g++/optim_slc3_323/stt.geant4-05-00-ref-03+tags1811/$ENV{'G4SYSTEM'}/$TestName";


#    $filename = "$ENV{'G4WORKDIR'}/stt.geant4-05-00-ref-03+tags1434/Linux-gO2/$TestName";
#    $filename = "/afs/cern.ch/user/s/stesting/stt/dev2/testtools/geant4/tests/tools/bin/icc/optim/Linux-icc/$TestName";

#    $filename = "$ENV{'G4WORKDIR'}/stt.geant4-05-00-cand-00+tags1024/$ENV{'G4SYSTEM'}/$TestName";
#    $filename = "$ENV{'G4WORKDIR'}/stt.geant4-05-00-cand-00+tags1035/$ENV{'G4SYSTEM'}/$TestName";
#    $filename = "$ENV{'G4WORKDIR'}/stt.geant4-04-01-ref-04+tags1019/$ENV{'G4SYSTEM'}/$TestName";

    $filename2 = "$SttDir2/$TestName";

    $errsize = (stat($filename))[7];

print "filename = $filename\n";
print "filename2 = $filename2\n";

    open(ERRFILE,"$filename") ;#|| die "Failed to open read $filename $!";
    open(ERRFILE2,"$filename2") ;#|| die "Failed to open read $filename $!";

    
    while(<ERRFILE>) {
	if ($ENV{'G4SYSTEM'} eq 'Linux-g++' || $ENV{'G4SYSTEM'} eq 'Linux-egcs') {
	    
#	    print "OK\n";
# It's very dependent on CPU load! Need to use user+system.
#	    if (/(\d+)\:(\d+)\.(\d+)elapsed/) {
	    if (/(\d+)\.(\d+)user/) {
		
		$usertime = "$1.$2";
		
		print "usertime=$usertime\n";
#		print "usertime-1=",$usertime-1,"\n";
	    } else {
#		$usertime=0;
	    }

	    if (/(\d+)\.(\d+)system/) {
		
		$systemtime = "$1.$2";
		
		print "systemtime=$systemtime\n";
#		print "systemtime-1=",$systemtime-1,"\n";
	    } else {
#		$systemtime=0;
	    }
	    
	    $deltatime1 = $usertime + $systemtime;
#	    print "$deltatime1\n";
	    if ( $deltatime1 ne 0) {
###		print "$TestName: deltatime1 = ",$deltatime1,"\n";
	    }

	} else {
	    if (/^real(\s+)(\d+)\.(\d+)/) {
		
		$deltatime1 = "$2.$3";
		
#		print "deltatime1=$deltatime1\n";
	    }
	}
    }

close(ERRFILE);

    while(<ERRFILE2>) {
	if ($ENV{'G4SYSTEM'} eq 'Linux-g++' || $ENV{'G4SYSTEM'} eq 'Linux-egcs') {
	    
#	    print "OK\n";
# It's very dependent on CPU load! Need to use user+system.
#	    if (/(\d+)\:(\d+)\.(\d+)elapsed/) {
	    if (/(\d+)\.(\d+)user/) {
		
		$usertime = "$1.$2";
		
#		print "usertime=$usertime\n";
#		print "usertime-1=",$usertime-1,"\n";
	    } else {
#		$usertime=0;
	    }

	    if (/(\d+)\.(\d+)system/) {
		
		$systemtime = "$1.$2";
		
#		print "systemtime=$systemtime\n";
#		print "systemtime-1=",$systemtime-1,"\n";
	    } else {
#		$systemtime=0;
	    }
	    
	    $deltatime2 = $usertime + $systemtime;
#	    print "$deltatime2\n";
	    if ( $deltatime2 ne 0) {
###		print "$TestName: deltatime2 = ",$deltatime2,"\n";
	    }

	} else {
	    if (/^real(\s+)(\d+)\.(\d+)/) {
		
		$deltatime2 = "$2.$3";
		
#		print "deltatime2=$deltatime2\n";
	    }
	}
#		print "$TestName: deltatime2 = ",$deltatime2,"\n";
    }
    
    close(ERRFILE2);



no integer;
$i = 0.29;
$j = 0.35;
$testNumber = $i/$j;

#print "testNumber = $testNumber\n";
  
#$ratio = $deltatime1/$deltatime2;

		print "$TestName: deltatime1 = ",$deltatime1,"\n";
		print "$TestName: deltatime2 = ",$deltatime2,"\n";

$TestName =~ s/.err//;
$Ratio = ($deltatime1*100)/$deltatime2;

use integer;
$IntRatio = $Ratio*1;
print "Ratio($TestName):  \t\t",$IntRatio,"\n";

}
    
