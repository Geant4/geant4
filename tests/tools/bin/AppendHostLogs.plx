#!/usr/local/bin/perl
#
#  AppendHostLogs.plx
#
$JunkHeap="/afs/cern.ch/sw/geant4/stt/testlogs/hostlogs";

@hosts=("dxplus","hpplus","pcgeant","refsol7","sungeant");

chomp($pwd=`pwd`);

print "Append and discard host log files\nfrom: $pwd\nto:   $JunkHeap\n";
opendir(HERE,"$pwd");
@files=readdir(HERE);
closedir(HERE);
foreach $file (@files) {
    if ( $file=~ m/^(\w+)\.(dev1|dev2|prod)\.(deb|opt).*\.log$/ ) {
        print "$file --> the JunkHeap\n";
# check $1 against @hosts
        open(FL,"$file") || die "Failed to open read $file\n";
        open(JUNK,">>$JunkHeap/$file") || die "Failed to open append $JunkHeap/$file $!";
        print JUNK "###\n### appending from $pwd\n";
        while ($line=<FL>) {
            chomp($line);
            print JUNK "$line\n";
        }
        close(JUNK);
        close(FL);
        unlink($file);
    }
}
