#!/usr/local/bin/perl
#
#  Hasami.plx
#
$NextShot="bonsai000.sdb";
if (defined($ARGV[0])) {$NextShot="bonsai$ARGV[0].sdb"}
print "Prepare $NextShot to drive Geant4 System Tests\n";
unlink("bonsai.sdb");
#
#  Extract from the database, needs perl version on pcgeant2 which requires
#  files in /afs/cern.ch/user/s/stesting/webtools/bonsai (or a change of path).

system("rsh pcgeant2 \"( cd /afs/cern.ch/user/s/stesting/webtools/bonsai; /home/sadilov/webtools/local/bin/perl -w makesdb.pl )\" ");
system("cp -p  /afs/cern.ch/user/s/stesting/webtools/bonsai/stt-dev.sdb bonsai.sdb");
unless(-e "bonsai.sdb") { 
    print "Unable to retrieve stt-dev.sdb from remote sql query on pcgeant2 \n"; 
    die;
}
$Age=(-M "bonsai.sdb"); # this might be fun
if ( $Age > 0 ) {
   print "Suspect a problem with sql query on pcgeant2,\n";
   print "  possible a clock synchronisation error\n";
   print "  will carry on regardless\n";
}
# Read the freshly generated bonsai.sdb
open(BONSAI,"bonsai.sdb") || die "Failed to open read bonsai.sdb $!";
# extract the reference tag (e.g. geant4-01-01-ref-05 )
$line = <BONSAI>;
if ( $line =~ m/^\.\s+(geant4-\d+-\d+-ref-\d+)/ ) {
    print "Last reference tag was $1\n" if ($verbose);
    $reftag=$1;
} else { 
    print "bonsai.sdb does not start with reference global tag\n";
    exit(1);
}
$reftagsdb="$reftag.sdb";
open(TAGGED,"$reftagsdb") || die "Failed to open read $reftagsdb $!";
while($done=<TAGGED>) {
    chomp($done);
    ($dir,$tag)=split(' ',$done);
    $dir=substr($dir,2);
    $InRefTag{$tag}=$reftag;
}
close(TAGGED);

while ($line = <BONSAI>) {
    chomp($line);
    ($dir,$tag)=split(' ',$line);
    $dir=substr($dir,2);            # neater without the leading ./
    printf ("%-25s ",$tag)if ($verbose);
    if ( defined($InRefTag{$tag})) {
        print " Done $InRefTag{$tag}" if ($verbose);
#       next;                       # we have accepted tags also superceeded in reftag.
    } 
    if ( defined($InThisRequest{$dir})) {
        print " Superceeds $InThisRequest{$dir}" if ($verbose);
        $xtag=$InThisRequest{$dir};
        $InThisTag{$xtag}=$dir;
    }
    $InThisRequest{$dir}=$tag;
    print "\n"if ($verbose);
}
close(BONSAI);

open(NEXT,">$NextShot") || die "Failed to open write next bonsai.sdb $NextShot $!";
open(BONSAI,"bonsai.sdb") || die "Failed to open read bonsai.sdb on second pass $!";
while ($line = <BONSAI>) {
    chomp($line);
    ($dir,$tag)=split(' ',$line);
    $dir=substr($dir,2);            # neater without the leading ./
    printf ("%-25s ",$tag);
    if ( defined($InRefTag{$tag})) { 
        print "     Done $InRefTag{$tag}\n";
        next; 
    } 
    if ( defined($InThisTag{$tag})) { 
        print "     superceeded\n";
        next;
    }
    print " $line\n";
    print NEXT "$line\n";
}
close(NEXT);
close(BONSAI);

print "\n\n  Enjoy and check your file $NextShot...\n\n";
open(NEXT,"$NextShot") || die "Failed to open read next bonsai.sdb $NextShot \n Just after I wrote it boys that is sick $!";
while ( <NEXT> ) { print }
close(NEXT);
