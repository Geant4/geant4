#!/usr/local/bin/perl
#
# /afs/cern.ch/sw/geant4/stt/dev1/testtools/geant4/tests/tools/bin
#
# Increment.plx
#
# Extract your new bonsai.sdb and other files OnTest
# update the source code (incrmentally if you wish).
#
# This script will then leave the tmp lib and bin directories in
# place and prepare a new stt(results) directory.
#
# Then run the tests incrementally and gmake will only do what is
# needed (maye).
#

$ActiveExamination="doit";
undef($ActiveExamination);
# ActiveExamination not coded below

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
#   next if ( "$host" eq "pcg4speed" );
    next if ( "$host" eq "pcgeant" );
    $G4SttWorkdir="/afs/cern.ch/sw/geant4/stt/$DevDir/$platform/$option";
    print "\nHost:      $host   \nWorkDir:      $G4SttWorkdir\n";
    $oldlink=`ssh $host "/bin/ls -l $G4SttWorkdir/stt "`;
    print "oldlink $oldlink\n";
    if ( $oldlink =~ m/->\s+(stt\.\S+)\s*/ ) {
        print "Previous stt directory was $1\n";
    } else {
        print "Previous stt directory not linked to stt\n";
    }
    $oldstat=`ssh $host "/bin/ls -l $G4SttWorkdir/*.stat "`;
    print "oldstat $oldstat\n";
    if ( $oldstat =~ m#/(\w+)\.stat\s*# ) {
        print "Previous status was $1\n";
        if ( "done" ne "$1" ) { 
            print "CHECK status $1 $G4SttWorkdir\n";
            system("ssh $host cat  $G4SttWorkdir/*.stat");
            print "NEXT Please\n";
            next;
        }
    } else {
        print "Previous status was not defined\n";
    }
    $checknew=`ssh $host "(/bin/ls -l $G4SttWorkdir/$SttTag >stdoutfile ) >& stderrfile ; cat stdoutfile stderrfile "`;
    unless ( $checknew =~ m/o\s+such\s+file\s+or/ ) {
        print "Requested new/incremented stt dir already exists\n";
        print "Skip      $host $G4SttWorkdir $SttTag\n\"$checknew\"\n";
        next;
    }
# 
# 
#     print "Trial csh redirection \"$checknew\"\n";
# 
# #    $checknew=`ssh $host "/bin/ls -l $G4SttWorkdir/$SttTag "`;
# #    $checklen=length($checknew);
# # maybe I need a csh manual, 2>&1 results ina null $checknew so for now
# # so maybe it runs in csh and not sh ?
# # just look at the length of the reply string on stdout rather than examine stderr
# #   unless ( $checknew =~ m/No\s+such/ || $checknew =~ m/not\s+found/ ) {
#     if ( $checklen > 6 ) {
#         print "Requested new/incremented stt dir already exists\n";
#         print "Skip      $host $G4SttWorkdir $SttTag\n\"$checknew\"\n";
#         next;
#     } else {
#         print "OK new directory does not already exist\n";
#     }
    print "Prepare $host $G4SttWorkdir $SttTag\n";

    system("ssh $host mkdir    $G4SttWorkdir/$SttTag"); # make (faked incremented) stt dir
    system("ssh $host mkdir    $G4SttWorkdir/$SttTag/$platform"); # required in incremental mode 
    system("ssh $host touch    $G4SttWorkdir/$SttTag/$platform/gmake.log"); # maybe echo later
    system("ssh $host /bin/rm  $G4SttWorkdir/*.stat");  # remove stat file (might lock)
    system("ssh $host /bin/rm  $G4SttWorkdir/stt");     # remove stt  link
    system("ssh $host \"(cd $G4SttWorkdir; /bin/ln -s  $SttTag stt)\" ");
}
close(TEMPLATE);
