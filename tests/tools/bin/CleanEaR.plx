#!/usr/local/bin/perl
#
#  CleanEaR.plx
#
open(CONFIG,"OnTest") || die "Failed to open OnTest configuration file $! ";
($DevDir,$Tag)=split(' ',<CONFIG>);
close(CONFIG);
$SttTag="stt." . "$Tag";
$TagNum="$Tag";
$TagNum=~s/^.*\+//;
print "Cleaning test $TagNum, working in dir \"$DevDir\"\nwith test-set tag \"$SttTag\"\n";
$TestTop="/afs/cern.ch/sw/geant4/stt/$DevDir";

opendir(TT,"$TestTop") || die "Failed to opendir TestTop  $TestTop $!";
@Platforms=grep(m/^[A-Z]/,readdir(TT));
closedir(TT);
foreach $Platform (@Platforms) {
   next unless (-d "$TestTop/$Platform");
#  print "$Platform\n";
   $PDir="$TestTop/$Platform";
   opendir(PLATFORM,"$PDir")  || die "Failed to opendir $PDir $! ";
   @Options=grep(!m/^\.\.?$/,readdir(PLATFORM));
   closedir(PLATFORM);
   foreach $Option (@Options) {
      $WorkDir="$PDir/$Option";
      opendir(WD,"$WorkDir") || { next };

# quick fix (Edit and Run)....

      print "Remove stt $SttTag directory and files\n";

# most likely we want the next line to execute.

#     system("rm  -rf  $WorkDir/$SttTag");


#     system("rm $PDir/$Option/tmp/$Platform/test18/exe/*.d");


      $LinkTo=`ls -l $WorkDir/stt`;
      chomp($LinkTo);
      unless ( $LinkTo =~m/.*-> (.*)$/ ) {next;}
      next unless ( "$1" eq "$SttTag");

      print "\nWe have tagset $1 in workdir\n    $WorkDir\n";
      @StatusFiles=grep(m/\.stat$/,readdir(WD));
      closedir(WD);
      foreach $StatusFile (@StatusFiles) {
          $StatusLine=`cat   $WorkDir/$StatusFile`;
          chomp($StatusLine);
          print "Remove $StatusFile (contents $StatusLine)\n";
          unlink("$WorkDir/$StatusFile");
      }
      print "Remove stt softlink\n";
      unlink("$WorkDir/stt");

  
    }
}
                      

