#!/usr/local/bin/perl
#
# /afs/cern.ch/sw/geant4/stt/dev1/testtools/geant4/tests/tools/bin
#
# Increment.plx
#
# update the source code then perform an incremental update of the 
# libraries and executables but write a new stt(results) directory
# then run the tests.

$ActiveExamination="doit";

open(CONFIG,"OnTest") || die "Failed to open OnTest configuration file $! ";
($DevDir,$Tag)=split(' ',<CONFIG>);
close(CONFIG);
print "\nStream:     \"$DevDir\"\n";
print "TagSet:     \"$Tag\"\n";
$SttTag="stt." . "$Tag";
$TestTop="/afs/cern.ch/sw/geant4/stt/$DevDir";
$TestLogDir="/afs/cern.ch/sw/geant4/stt/$DevDir/testtools/geant4/tests/tools/bin";

opendir(TL,"$TestLogDir") || die "Failed to opendir TestLog  $TestLogDir $!";
# @testlogs=grep(m/^\w+\.dev\d\.\w+.*\.log/,readdir(TL));
closedir(TL);

open(TEMPLATE,"stthosts.data") || die "Failed to open read stthosts.data  $!";
while ( <TEMPLATE> ) {
    chomp;
    if ( /^#/) { next}
    if ( /^%/) { next}
    ($rundir,$host,$platform,$option,@fiveargs)=split(' ');
    next if ( "$rundir" ne "$DevDir");
    $G4SttWorkdir="/afs/cern.ch/sw/geant4/stt/$DevDir/$platform/$option";
    print "\nHost:      $host   \nWorkDir:      $G4SttWorkdir\n";
    $oldlink=`rsh $host "/bin/ls -l $G4SttWorkdir/stt"`;
    if ( $oldlink =~ m/->\s+(stt\.\S+)\s*/ ) {
        print "Previous stt directory was $1\n";
    } else {
        print "Previous stt directory not linked to stt\n";
    }
    $oldstat=`rsh $host "/bin/ls -l $G4SttWorkdir/*.stat"`;
    if ( $oldstat =~ m#/(\w+)\.stat\s*# ) {
        print "Previous status was $1\n";
        if ( "done" ne "$1" ) { 
            print "CHECK status $1 $G4SttWorkdir\n";
            system("rsh $host cat  $G4SttWorkdir/*.stat");
        }
    } else {
        print "Previous status was not defined\n";
    }
    $checknew=`rsh $host "/bin/ls -l $G4SttWorkdir/$SttTag "`;
    $checklen=length($checknew);
# maybe I need a csh manual, 2>&1 results ina null $checknew so for now
# so maybe it runs in csh and not sh ?
# just look at the length of the reply string on stdout rather than examine stderr
#   unless ( $checknew =~ m/No\s+such/ || $checknew =~ m/not\s+found/ ) {
    if ( $checklen > 6 ) {
        print "Requested new/incremented stt dir already exists\n";
        print "Skip      $host $G4SttWorkdir $SttTag\n\"$checknew\"\n";
        next;
    }
    print "Prepare $host $G4SttWorkdir $SttTag\n";
    system("rsh $host mkdir    $G4SttWorkdir/$SttTag"); # make (faked incremented) stt dir
    system("rsh $host mkdir    $G4SttWorkdir/$SttTag/$platform"); # required in incremental mode 
    system("rsh $host touch    $G4SttWorkdir/$SttTag/$platform/gmake.log"); # maybe echo later
    system("rsh $host /bin/rm  $G4SttWorkdir/*.stat");  # remove stat file (might lock)
    system("rsh $host /bin/rm  $G4SttWorkdir/stt");     # remove stt  link
    system("rsh $host \"(cd $G4SttWorkdir; /bin/ln -s  $SttTag stt)\" ");
}
close(TEMPLATE);
