#!/usr/local/bin/perl


$SttDir = "$ENV{'G4WORKDIR'}/stt/$ENV{'G4SYSTEM'}/";

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
    $errsize = (stat($filename))[7];

    open(ERRFILE,"$filename") || die "Failed to open read $filename $!";
    
    while(<ERRFILE>) {
	if ($ENV{'G4SYSTEM'} eq 'Linux-g++' || $ENV{'G4SYSTEM'} eq 'Linux-egcs') {
# It's very dependent on CPU load! Need to use user+system.
#	    if (/(\d+)\:(\d+)\.(\d+)elapsed/) {
	    if (/(\d+)\.(\d+)user/) {
		
		$sec = $1*60+$2;
		$deltatime = "$sec.$3";
		
		print "deltatime=$deltatime\n";
	    }
	} else {
	    if (/^real(\s+)(\d+)\.(\d+)/) {
		
		$deltatime = "$2.$3";
		
#		print "deltatime=$deltatime\n";
	    }
	}
    }
    
    close(ERRFILE);

print "$TestName: deltatime=$deltatime\n";

}
    
