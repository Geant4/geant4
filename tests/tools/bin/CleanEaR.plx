#!/usr/local/bin/perl
#
# /afs/cern.ch/sw/geant4/stt/dev2
#
#  CleanEar.plx
#
# /afs/cern.ch/sw/geant4/stt/dev1/Linux-g++/optim/tmp/Linux-g++/test18/exe
#
open(CONFIG,"OnTest") || die "Failed to open OnTest configuration file $! ";
($DevDir,$Tag)=split(' ',<CONFIG>);
close(CONFIG);
print "Working in dir \"$DevDir\" with test-set tag \"$Tag\"\n";
$SttTag="stt." . "$Tag";
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
      opendir(WD,"$WorkDir") || { print  "Failed to opendir WorkDir $WorkDir $! \n"} ;
      opendir(WD,"$WorkDir") || { next };

# quick fix (Edit and Run)....

      system("ls -l $WorkDir/stt");
#     system("rm  $WorkDir/stt");
#     system("ls -l $WorkDir/*stat");
      system("cat   $WorkDir/*stat");
#     system("rm    $WorkDir/*stat");

      system("rm $PDir/$Option/tmp/$Platform/test18/exe/*.d");
  
      $TestDir="$PDir/$Option/$SttTag/$Platform";
   #  opendir(TD,"$TestDir") || { print  "Failed to opendir TestDir $TestDir $! \n"} ;
      opendir(TD,"$TestDir") || { next };
    }
}
                      

