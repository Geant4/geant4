#!/usr/local/bin/perl
#
# /afs/cern.ch/sw/geant4/stt/dev1/testtools/geant4/tests/tools/bin
# ExtractLastLog.plx
#

open(CONFIG,"OnTest") || die "Failed to open OnTest configuration file $! ";
($DevDir,$Tag)=split(' ',<CONFIG>);
close(CONFIG);
print "\nWorking in dir \"$DevDir\" with test-set tag \"$Tag\"\n";
$SttTag="stt." . "$Tag";
$TestTop="/afs/cern.ch/sw/geant4/stt/$DevDir";
$TestLogDir="/afs/cern.ch/sw/geant4/stt/$DevDir/testtools/geant4/tests/tools/bin";

opendir(TL,"$TestLogDir") || die "Failed to opendir TestLog  $TestLogDir $!";
@testlogs=grep(m/^\w+\.dev\d\.\w+.*\.log/,readdir(TL));
closedir(TL);
foreach $testlog (@testlogs) {
    next unless ((-M "$TestLogDir/$testlog") < 60 );
    print "Test Run Log \"$testlog\"\n";
    $Machine="Machine";
    $Option="CompilerOpts";
#  dxplus04.dev1.opt_NONISO.log
    if ( $testlog =~ m/(\w+)\.\w+\.(.*)\.log/ ) {
        $Machine=$1;
        $Option=$2;
    }
    $lines=0;
    $title=0;
    $copy=0;
    open(TLC,"$TestLogDir/$testlog") || die "Failed to open read $TestLogDir/$testlog $!";
    while ($line = <TLC> ) {
        $lines++;
        if ( $line =~ /^_____/ ) {$copy=0;
            # next; Problem with metacharacter in tag - use index on next line for now
            $newtestline=$line;
            $line = <TLC>;
            $copy = (index($line,$Tag) > 0);
        }
#       if ( $line =~ /\$Tag/) {$copy=1;print "MATCH TAG <<<<<<<<<<<<<<\n";}
        if ( $copy ) {
            if ( 0 == $title ) { print "\n\nTag:  $Tag\nTest:  $testlog\n\n$newtestline";$title=2;}
            $lines++;
            chomp($line);
            print "$Machine $Option $lines    $line\n";
        }
    }
    close(TLC);
    print "\n";
}
exit();
