#!/usr/local/bin/perl
#
# ReviewUpdate.plx
#
open(CONFIG,"OnTest") || die "Failed to open OnTest configuration file $! ";
($DevDir,$Tag)=split(' ',<CONFIG>);
close(CONFIG);
print "Working in dir \"$DevDir\" with test-set tag \"$Tag\"\n";
if ( $Tag =~ m/tags(\d+)/ ) {
   print "Review update log for tagset $1\n";
   $UpdateLog="update$1.log.1";
   $UpdateLog="update$1.log";
} else {
   print "Failed to extract tagset number from OnTest configuration file\n";
   exit(1);
}
$CVSModifieds=0; 
$CVSConflicts=0;
$LogLevel=1;
$logfilesize=(-s $UpdateLog);
print "$UpdateLog is $logfilesize characters\n";
open(ULOG,"$UpdateLog") || die "Failed to open (read) $UpdateLog $!";
while ($line = <ULOG> ) {
    $lines++;
    chomp($line);
    if ( $line =~ /^cvs update / ) {
        $CVSCommand = $line;
        $CVSCommands++;
        print "Command  $line\n" if ( $LogLevel > 1);
        next;
    }
    if ( $line =~ /^U / ) {
        $CVSUpdated = $line;
        $CVSUpdated =~ s/^U //;
        $CVSFileUpdates++;
        next;
    }
    if ( $line =~ /^cvs update: Updating / ) {
        $CVSUpdating = $line;
        $CVSUpdating =~ s/.*Updating\s+//;
        next;
    }
    if ( $line =~ /^cvs update: (.*) is no longer in the repository/ ) {
        $CVSNotInRepository=$1;
        $CVSFileNotInRepository++;
        next;
    }
    if ( $line =~ /^\?/ ) {
        $Questionable++;
        if ( $line =~ /rndm$/ ) { $Passrndm++; next}
        if ( $line =~ /ir.out$/ ) { $Passirout++; next}
        if ( $line =~ /\.log$/ ) { $Passdotlog++; next}
        if ( $line =~ m#tests/tools/bin# ) { $Passtools++; next}
        print "Updating $CVSUpdating      Questionable file $line\n";
        next;
    }
# cvs update: warning: source/run/src/G4RunManager.cc was lost
    if ( $line =~ /cvs update: warning:\s+(.*)\s+was lost/ ) {
        $CVSLost++;
        if ( $1 =~ /G4RunManager.cc/ ) { $CVSLostOK++; next}
        print "Lost:   $line\n" if ( $LogLevel > 1);
        next;
    }
    if ( $line =~ /^C / ) {
        $CVSConflict = $line;
        $CVSConflict =~ s/^C //;
        $CVSConflicts++;
        print "$line\n" if ( $LogLevel > 1);
        print "Conflict: $CVSConflict\n" if ( $LogLevel > 0);
        push(@Conflict,$CVSConflict);
        next;
    }
    if ( $line =~ /^M / ) {
        $CVSModified = $line;
        $CVSModified =~ s/^C //;
        $CVSModifieds++;
        print "$line\n"  if ( $LogLevel > 1);
        print "Modified: $CVSModified\n" if ( $LogLevel > 0);
        push(@Modified,$CVSModified);
        next;
    }


    print "$line\n";
    $displayed++;
}
close(ULOG);

$QOk=$Passrndm+$Passirout+$Passtools+$Passdotlog;
print "\n\n";
print "Review update log for tagset $1\n";
print "Displayed:           $displayed lines of $lines\n";
print "Commands:            $CVSCommands\n";
print "File Updates:        $CVSFileUpdates\n";
print "Questionable:        $QOk of $Questionable ignored\n";
print "Not in Repository:   $CVSFileNotInRepository\n";
print "Was lost:            $CVSLostOK normal (G4RunManager.cc) of $CVSLost\n";
print "Modified:            $CVSModifieds  investigate if not 0\n";
print "Conflicts:           $CVSConflicts  investigate if not 0\n";

foreach (@Conflict) { print "$_\n"; }


