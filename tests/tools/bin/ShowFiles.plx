#!/usr/local/bin/perl
#
# /afs/cern.ch/sw/geant4/stt/dev2
#
#  ShowFiles.plx
#  The comparison part of this script may eventually be part of the
#  testing but much can be done withing the testing script and reported
#  during execution.
#
$DevDir="dev2";
$SttTag="stt.geant4-01-01-ref-03-tags40";
$DevDir="dev1";
$SttTag="stt.geant4-01-01-ref-04+group1a";
#
open(CONFIG,"OnTest") || die "Failed to open OnTest configuration file $! ";
($DevDir,$Tag)=split(' ',<CONFIG>);
close(CONFIG);
print "Working in dir \"$DevDir\" with test-set tag \"$Tag\"\n";
$SttTag="stt." . "$Tag";
$TestTop="/afs/cern.ch/sw/geant4/stt/$DevDir";

@filetypes=qw ( out err diff diff_err );

@tl = (
       [ qw ( test01    test02    test02.hadron    test05    test06     )],
       [ qw ( test07    test09    test10           test101   test102    )],
       [ qw ( test103   test104   test104.EMtest   test105   test106    )],
       [ qw ( test11    test12    test13           test14    test15     )],
       [ qw ( test16    test401   test402          test501   test502    )], 
       [ qw ( test503   test504   test505   test506)] 
      );

opendir(TT,"$TestTop") || die "Failed to opendir TestTop  $TestTop $!";
@Platforms=grep(m/^[A-Z]/,readdir(TT));
closedir(TT);
foreach $Platform (@Platforms) {
   print "$Platform\n";
#  next unless (-d "$TestTop/$Platform");
   $PDir="$TestTop/$Platform";
   next unless (-d "$PDir");
   opendir(PLATFORM,"$PDir")  || die "Failed to opendir $PDir $! ";
   @Options=grep(!m/^\.\.?$/,readdir(PLATFORM));
   closedir(PLATFORM);
#
# DEC, HP, AIX are in afs
# Linux and Sun are linked/mounted locally.
#
   foreach $Option (@Options) {
      $TestDir="$PDir/$Option/$SttTag/$Platform";
      opendir(TD,"$TestDir") || { print  "Failed to opendir TestDir $TestDir $! "} ;
      opendir(TD,"$TestDir") || { next };

      $Heading="$Platform  $Option  $SttTag";
      $lenH=length($Heading);
      $lenS=int((80-$lenH)/2);
      $Space=' ';
      $Under='=';
      $H1="$Space" x $lenS; $H1.="$Heading";
      $H2="$Space" x $lenS; $H2.="$Under" x $lenH;
      print "\n\n$H1\n";
      print "$H2\n";
# sort of debug
      $text1=`ls -l $PDir/$Option/*stat`;
      print "$text1";
      $text2=$text1;
      $text2=~s#^.*/afs#/afs#;
      print "$text2";
      $text3=`cat   $PDir/$Option/*stat`;
      print "$text3";
# sort of debug

      undef(%Filesize);
      @tfiles=readdir(TD);
      closedir(TD);
      $playing="analyse";
      foreach $tfile (@tfiles) {
         $fsize=(-s "$TestDir/$tfile");
         if ( $tfile =~ m/(test\d{1,3}\.??\w*)\.(\w+)$/ ) {
             $Filesize{$1}{$2}=$fsize;
             if ( $playing ) {
                 $TestName=$1;
                 $FileType=$2;
#                print "$Platform  $Option  $SttTag $TestName $FileType $fsize \n";
                 if ( "err" eq "$FileType" ) {
                     open(ANAERR,"$TestDir/$tfile");
                     while ($line=<ANAERR>) {
                         next if ($line =~m/^real\s+\d/);
                         next if ($line =~m/^user\s+\d/);
                         next if ($line =~m/^sys\s+\d/);
                         next if ($line =~m/user.*system.*elapsed.*%CPU/);
                         next if ($line =~m/nputs.*outputs.*pagefaults.*swaps/);
                         next if ($line =~m/^\s+$/);
                         chomp($line);
                         $report_err=sprintf("%-16s %-9s %-14s \"%s\".",$TestName,$Platform,$Option,$line);
#                        print "ERROR: $report_err\n";
                         push(@ErrorReports,$report_err);
                     }
                     close(ANAERR);
                 }
                 if ( "diff" eq "$FileType" ) {
                     $maxshow=10; $shown=0; $differs=0;
                     open(ANADIFF,"$TestDir/$tfile");
                     while ($line=<ANADIFF>) {
                         $differs++;
                         next if ($line =~m/^\d+/);                        # diff generates ed linemode command
                         next if ($line =~m/^\s+$/);                       # blank 
                         next if ($line =~m/^---$/);                       # diff marks end of compared sections
                         next if ($line =~m/Geant4 version \$Name:/);      # CVS Tags
#                        next if ($line =~m/\$Id: ShowFiles.plx,v 1.2 2000-06-05 14:54:42 stesting Exp $/);               # CVS Tags
                         next if ($line =~m/\$Id:/);               # CVS Tags
                         next if ($line =~m/\s+\(\d{2}-[A-Z][a-z]{2}-\d{2,4}\)/);
                         next if ($line =~ m/User=.*Real=.*Sys=/);
                         next if ($line =~ m/Material:.* 6317232936888469/);     # ISO printing
                         next if ($line =~ m/\s+\(pointer\)\s+/);          # pointers will surely be different
                         next if ($line =~ m/\s+pointer\s+/);              # pointers will surely be different
                         next if ($line =~ m#/afs/cern.ch/sw/geant4/stt/dev\d+/src/geant4/#);   # tests in dev1 or dev2 
                         next if ($line =~ m/Current couple of seeds =/);  # seeds will eventually be different
                         $shown++;
                         next if ( $shown > $maxshow );
                         chomp($line);
                         print "$TestName $FileType  \"$line\".\n";
                     }
                     close(ANADIFF);
                 #   print "$TestName $FileType $shown significant differences\n";
                     if ($shown > 0 ) {
                         $report_diff=sprintf("%-16s %-9s %-14s %7d different lines (in diff file).",$TestName,$Platform,$Option,$shown);
                         push(@DiffReports,$report_diff);
                     }
                 }
             }
         }
      }

      undef(%Total),undef(%Lineout);
      foreach $i ( 0 .. $#tl ) {
         print "\n\n      $Platform  $Option  $SttTag\n";
         print "\nTest    :";
         foreach $test ( @{$tl[$i]} ) {
            printf ("%15s",$test);
         }
         if ( $i == $#tl ) {printf ("%20s","Total");}
         foreach $step ("err","diff_err","diff","out") {
              printf ("\n%8s:",$step);
          
              foreach $test ( @{$tl[$i]} ) {
                  printf ("%15d", $Filesize{$test}{$step} );
                  $Total{$step}+=$Filesize{$test}{$step};
                  $Lineout{$test}+=$Filesize{$test}{$step};
              }
              if ( $i == $#tl ) {printf ("%20d",$Total{$step});}
         }
         print "\n";
      }
      foreach $test (keys(%Lineout)) {
          next unless ( $Lineout{$test} <= 0 );
          $report_notrun=sprintf("%-16s %-9s %-14s No output files from this test.",$test,$Platform,$Option);
#         print "NOTRUN: $report_notrun\n";
          push(@NotRunReports,$report_notrun);
      }
   }
}
if (defined(@NotRunReports)) {
   print "\n\n\n                 Tests NOT RUN in $SttTag\n\n";
   foreach $line (sort(@NotRunReports)) {
       print "$line\n";
   }
}
if (defined(@ErrorReports)) {
   $sameagain="no";
   print "\n\n\n                 ERRORS seen in $SttTag\n\n";
   foreach $line (sort(@ErrorReports)) {
       $same=substr($line,0,16);
       print "\n" if ("$same" ne "$sameagain" );
       print "$line\n";
       $sameagain=$same;
   }
}
if (defined(@DiffReports)) {
   print "\n\n\n                 Different lines from diff in $SttTag\n";
   $sameagain="no";
   foreach $line (sort(@DiffReports)) {
       $same=substr($line,0,16);
       print "\n" if ("$same" ne "$sameagain" );
       print "$line\n";
       $sameagain=$same;
   }
}
