#!/usr/local/bin/perl


$SttDir = "$ENV{'G4WORKDIR'}/stt.geant4-04-01-ref-00+tags949/$ENV{'G4SYSTEM'}/";
$SttDir2 = "$ENV{'G4WORKDIR'}/stt.geant4-04-01-ref-00+tags935/$ENV{'G4SYSTEM'}/";


#$SttDir2 = $SttDir;
#$SttDir2 =~ s|prod|dev1|;

print "SttDir  = $SttDir\n";
print "SttDir2 = $SttDir2\n";

my $deltatime1 ;
my $deltatime2 ;

opendir(TT,"$SttDir") || die "Failed to opendir  $SttDir $!";
@Tests = grep(m/\.err$/,readdir(TT));
closedir(TT);

foreach $Test (@Tests) {
#

$TestName = $Test;

print "TestName = $TestName\n";

#
# Get the size of $TestName.err 
#
    $filename = "$ENV{'G4WORKDIR'}/stt/$ENV{'G4SYSTEM'}/$TestName";
    $filename2 = "$SttDir2/$TestName";

    $errsize = (stat($filename))[7];

print "filename = $filename\n";
print "filename2 = $filename2\n";

    open(ERRFILE,"$filename") || die "Failed to open read $filename $!";
    open(ERRFILE2,"$filename2") || die "Failed to open read $filename $!";

    
    while(<ERRFILE>) {
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
	    
	    $deltatime1 = $usertime + $systemtime;
#	    print "$deltatime1\n";
	    if ( $deltatime1 ne 0) {
		print "$TestName: deltatime1 = ",$deltatime1,"\n";
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
		print "$TestName: deltatime2 = ",$deltatime2,"\n";
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

print "testNumber = $testNumber\n";
  
#$ratio = $deltatime1/$deltatime2;

		print "$TestName: deltatime1 = ",$deltatime1,"\n";
		print "$TestName: deltatime2 = ",$deltatime2,"\n";

$TestName =~ s/.err//;
$Ratio = ($deltatime1*100)/$deltatime2;

use integer;
$IntRatio = $Ratio*1;
print "Ratio($TestName): new/release = \t\t",$IntRatio,"\n";

}
    
