#!/usr/local/bin/perl
#
#  Hasami.plx
#
$TagsDir="/afs/cern.ch/sw/geant4/stt/testlogs/tagsets";
require("/afs/cern.ch/user/s/stesting/scripts/Services/G4SttMonitor.pl");


$devprod=&CurrentSttDir();

if (defined($ARGV[0])) {
    $tagset=$ARGV[0]
} else {
    $tagset=&NextTagSetNumber()
}

$NextShot="bonsai$tagset.sdb";
$NextLog="pruned$tagset.log";

print "Prepare Geant4 System Tests in $devprod for tagset $tagset\n";

unlink("bonsai.sdb");
&ExtractBonsaiSdb();

# Read the freshly generated bonsai.sdb
# Some operations are hard coded and others are from other files. 
# Overwrite geant4-01-00-cand-00 with geant4-02-00-cand-00

open(BONSAI,"bonsai.sdb") || die "Failed to open read bonsai.sdb $!";

# extract the reference tag (e.g. geant4-01-01-ref-05 )

$line = <BONSAI>;
if ( $line =~ m/geant4-01-00-cand-00/ ) {
    $line =~ s/geant4-01-00-cand-00/geant4-02-00-cand-00/;
    print "Fixup for database problem\n$line\n";
}
if ( $line =~ m/^\.\s+(geant4-\d+-\d+-(ref|cand)-\d+)/ ) {
    print "Last reference tag was $1\n" if ($verbose);
    $reftag=$1;
} else { 
    print "bonsai.sdb does not start with a recognised reference global tag\n";
    exit(1);
}

# read the sdb file with the tags included in the reference tag.

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
open(PLOG,">$NextLog") || die "Failed to open write pruning log $NextLog $!";
open(BONSAI,"bonsai.sdb") || die "Failed to open read bonsai.sdb on second pass $!";
while ($line = <BONSAI>) {
    chomp($line);
    if ( $line =~ m/geant4-01-00-cand-00/ ) {
        $line =~ s/geant4-01-00-cand-00/geant4-02-00-cand-00/;
        print "Fixup for database problem\n$line\n";
        print PLOG "Fixup for database problem\n$line\n";
    }
    ($dir,$tag)=split(' ',$line);
    $dir=substr($dir,2);            # neater without the leading ./
    printf PLOG ("%-25s ",$tag);
    if ( defined($InRefTag{$tag})) { 
        print PLOG "     Done $InRefTag{$tag}\n";
        next; 
    } 
    if ( defined($InThisTag{$tag})) { 
        print PLOG "     superceeded\n";
        next;
    }
    print PLOG " $line\n";
    print NEXT "$line\n";
}
close(NEXT);
close(BONSAI);
close(PLOG);

$hostdata = &CreateHostsFile( $reftag , $tagset, $devprod );

open(ONTEST,">OnTest") || die "Failed to open (write) OnTest $!";
print ONTEST "$devprod ${reftag}+tags${tagset}\n";
close(ONTEST);

print "\nUpdate your code with          $NextShot .....\n";
print "The cuttings are in            $NextLog\n";
print "The tests may be driven with   $hostdata\n";
print "OnTest we now have             $devprod ${reftag}+tags${tagset}\n";

exit();


sub ExtractBonsaiSdb {
    
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
}

sub NextTagSetNumber {
    my($Service);
    my($TagSetNumberFile);
    my($TagSetNumber);
    $TagSetNumberFile="$TagsDir/tagsetnumber";
    
    $Service="NextTagSetNumber";
    &G4SttMonitor("$Service",'lock',10)  && die "$Service Locked";
    open(TSN,"$TagSetNumberFile") || die "Failed to open (read) $TagSetNumberFile $!";
    $TagSetNumber=<TSN>;
    chomp($TagSetNumber);
    close(TSN);
    $TagSetNumber++;
    open(TSN,">$TagSetNumberFile") || die "Failed to open (write) $TagSetNumberFile $!";
    print TSN "$TagSetNumber\n";
    close(TSN);
    &G4SttMonitor("$Service",'free',10)  && die "$Service was not Locked at end of Servicing";
    return $TagSetNumber;
}
sub CurrentSttDir {
    my($here);
    my($devprod);
    $here=`pwd`;
    if ( $here =~ m#/((dev1|dev2|prod))/test# ) {
        $devprod=$1;
    } else {
        die "These scripts are bound to the dev1 or dev2 or prod structures";
    }
    return $devprod; # I hate returning a value from the middle of the code
}
sub CreateHostsFile {
#   arguments $reftag , $tagset, $devprod 
    my($sttdir,$hostdata,$host,$platform,$option,@fiveargs);
    $sttdir="$_[0]+tags$_[1]";
    $hostdata="hosts$_[1].data";
    open(HOSTFILE,">$hostdata") || die "Failed to open (write) $hostdata  $!";
    open(TEMPLATE,"stthosts.data") || die "Failed to open read stthosts.data  $!";
    while ( <TEMPLATE> ) {
        chomp;
        if ( /^#/) { print HOSTFILE "$_\n"; next}
        if ( /^%/) { print HOSTFILE "$_\n"; next}
        ($host,$platform,$option,@fiveargs)=split(' ');
        printf HOSTFILE ("%-8s %-4s  %-10s %s  %s %s %s %s %s\n",$host,$_[2],$option,$sttdir,@fiveargs);
    }
    close(HOSTFILE);
    close(TEMPLATE);
    return $hostdata;
}
