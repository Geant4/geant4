#!/usr/local/bin/perl
#
# /afs/cern.ch/sw/geant4/stt/dev1/testtools/geant4/tests/tools/bin
# ExtractLastLog.plx
#

$ActiveExamination="doit";

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
#   print "Test Run Log \"$testlog\"\n";
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
    undef($Step); undef($ResultsDir);
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
            if ( $ShowAll ) {print "$Machine $Option $lines    $line\n"};
            if ( $line =~ m/^STT:(\w+)\s+(\w+)/ ) {
                print "$Machine $Option $lines    $line\n";
                $Step=$1; $State=$2;
# Collect subsidiary information
                if ( $line =~ m/^STT:BUILD Started/ ) {
                     $nextline = <TLC>;
                     if ( $nextline =~ m#^(/afs/.*)\s+created# ) {
                         $ResultsDir=$1;
                         undef($syntaxerror);
                     } else {
                         print "Unexpected text in log file\n$nextline\ncreation of workdir expected\n";
                     }
                }
            }
            if ( $line =~ m/^Starting (test\d+) in/ ) { $Start{$1}=$line }
            if ( $line =~ /Finished (test\d+) in/ ) { $Finish{$1}=$line }
            if ( $line =~ /Disabled (test\d+) in/ ) { $Disabled{$1}=$line }
            if ( $line =~ m/syntax\s+error/ ) {
                $syntaxerror=$line;
                print "syntax error reported\n$line\n";
            }
        }
    }
    close(TLC);
    next unless(defined($Step));
    print "Reached $Step $State Results Directory is ...\n";
    print "$ResultsDir\n";
    print "\n";
    next if ( $Machine =~ /dxplus00/ );
    if ( $ActiveExamination && "$Step.$State" eq "BUILD.Started" ) {
#       @MakingProgress=`rsh -l stesting $Machine "ls -l $ResultsDir/gmake.log"`;
#       foreach $line (@MakingProgress) {
#           chomp($line);
#           print "Progress $line\n";
#       }
        @MakingProgress=`rsh -l stesting $Machine "/afs/cern.ch/sw/geant4/stt/dev1/testtools/geant4/tests/tools/bin/FileAge.plx $ResultsDir/gmake.log"`;
        foreach $line (@MakingProgress) {
            chomp($line);
            print "Gmake Age $line\n";
        }
    }
    if (defined($syntaxerror)) {
        print "$syntaxerror\n";
        print "SE-- $ResultsDir\n";
    }
    undef($pl);undef($sl);
    foreach $testnum (keys (%Start)) {
        $num=$testnum;
        $num=~s/test//;
        $TestState{$testnum}="S";
        if (defined($Finish{$testnum})) {$TestState{$testnum}.="F"}
        $pl.=sprintf("%4d",$num);
        $sl.=sprintf("%4s",$TestState{$testnum});
    }
    print "Test  $pl\n";
    print "State $sl\n";
    
}
exit();
