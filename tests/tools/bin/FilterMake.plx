#!/usr/local/bin/perl
#
# /afs/cern.ch/sw/geant4/stt/dev2
#
#  FilterMake.plx
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
      $TestDir="$PDir/$Option/$SttTag/$Platform";
#     opendir(TD,"$TestDir") || { print  "Failed to opendir TestDir $TDir $! "} ;
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

      undef(%Filesize);
      @tfiles=readdir(TD);
      closedir(TD);
      @tfiles=("gmake.log");
      foreach $tfile (@tfiles) {
         $fsize=(-s "$TestDir/$tfile");
         if ( $tfile =~ m/(gmake)\.(\w+)$/ ) {
             $Filesize{$1}{$2}=$fsize;
             $TestName=$1;
             $FileType=$2;
             print "$Platform  $Option  $SttTag $TestName $FileType $fsize \n";
             undef($playing );
              
             open(GMAKE,"$TestDir/gmake.log");
             while ($line=<GMAKE>) {
                 chomp($line);
                 if ( $line =~ m#Making dependency for file (\w+)\/(\w+)\.#) {
                     $Action="Dependency";
                     $Routine="$1/$2";
                     next;
                 }
                 if ( $line =~ m#Making dependency for file (\w+)\.#) {
                     $Action="Dependency";
                     $Routine="$1";
                     next;
                 }
                 if ( $line =~ m#Compiling (\w+)\.#) {
                     $Action="Compiling";
                     $Routine="$1";
                     next;
                 }
                 if ( $line =~ m#Creating/replacing object files in .*\/(\w+\/\w+\.a)#) {
                     $Action="CreateLib";
                     $Routine="$1";
                     next;
                 }
                 if ( $line =~ m#Creating/replacing object files in .*\/(\w+\.a)\s*#) {
                     $Action="CreateLib";
                     $Routine="$1";
                     next;
                 }
                 if ( $line =~ m#^CC # || $line =~ m#^aCC #) {
                     $Action="cxx";
                     $Routine="cxx command";
                     next;
                 }
                 if ( $line =~ m#^cxx # || $line =~ m#^g\+\+ #) {
                     $Action="cxx";
                     $Routine="cxx command";
                     next;
                 }
                 if ( $line =~ m#^f77 # ) {
                     $Action="f77";
                     $Routine="Fortran77";
                     next;
                 }
                 if ( $line =~ m#^cxx: Warning:#) {
                     $Action="C++ warning";
                     $Routine="cxx command";
                     next;
                 }
                 if ( $line =~ m#^Warning 652:#) {
                     $Action="C++ warning";
                     $Routine="cxx command";
                     next;
                 }
                 if ( $line =~ m#^cd \.\./tools/lib\; gmake#) {
                     next;
                 }
# C++ Warning from SUN-CC
                 if ( $line =~ m#^\"\w+\/(\w+\.cc)\"\,\s+line\s+(\d+)\:\s+(\w+)\:(.*\.)# ) {
#                    print "SUN C++ file $1 line $2 severity $3 text $4\n";
                     next;
                 }
                 if ( $line =~ m#^\"\w+\/(\w+\.hh)\"\,\s+line\s+(\d+)\:\s+(\w+)\:(.*\.)# ) {
#                    print "SUN Header file $1 line $2 severity $3 text $4\n";
                     next;
                 }

# We should read these lines from a file.

                 if ( $line =~ m#^WARNING: Making a library map of granular libraries.# ) { next; }
                 if ( $line =~ m#^         This is a list of libraries in order of use, and for# ) { next; }
                 if ( $line =~ m#^         each library a list of other libraries used.# ) { next; }
                 if ( $line =~ m#^         To do this it needs a complete set of dependency# ) { next; }
                 if ( $line =~ m#^         files, e.g., after gmake in the source/ directory.# ) { next; }
                 if ( $line =~ m#^Searching /afs/cern.ch/sw/geant4/stt/dev2/src/geant4/source# ) { next; }
                 if ( $line =~ m#^  for GNUmakefiles containing "name" and sorting...# ) { next; }
                 if ( $line =~ m#^Weeding out global level GNUmakefiles and non-libraries...# ) { next; }
                 if ( $line =~ m#^Making libname.map starter file...# ) { next; }
                 if ( $line =~ m#^Making libname.map...# ) { next; }
                 if ( $line =~ m#^  Reading library name map file...# ) { next; }
                 if ( $line =~ m#^  Reading dependency files...# ) { next; }
                 if ( $line =~ m#^  Checking for circular dependencies...# ) { next; }
                 if ( $line =~ m#^  Reordering according to dependencies...# ) { next; }
                 if ( $line =~ m#^  Writing new library map file...# ) { next; }

                 if ( $line =~ m#^real\s+\d+\.\d+# ) { next; }
                 if ( $line =~ m#^user\s+\d+\.\d+# ) { next; }
                 if ( $line =~ m#^sys\s+\d+\.\d+# ) { next; }
                 if ( $line =~ m#^real\s+\d+\:\d+# ) { next; }
                 if ( $line =~ m#^user\s+\d+\:\d+# ) { next; }
                 if ( $line =~ m#^sys\s+\d+\:\d+# ) { next; }



                 if ( $line =~ m#^gmake all#) {
                     next;
                 }
                 if ( $line =~ m#^gmake\[.*is up to date\.$#) {
                     next;
                 }
                 if ( $line =~ m#^ar: Warning: creating#) {next}
                 if ( $line =~ m#^ar: creating#) {next}            # hp aCC
                 if ( $line =~ m#^Using granular#) {next}
                 if ( $line =~ m#^Linking#) {next}

#                if ( "$Action" eq "C++ warning" ) {next}
                 if ( $line =~ m#^ild:# ) { $ild++; next}
                 if ( $line =~ m#^\"src/SdaiCONFIG_CONTROL_DESIGN.cc\", line (\d+): (\w+): (.*)$# ) {
                     $LineNum=$1; $Severity=$2; $Whinge=$3;
                     print " SDAI: $LineNum $Severity $Whinge\"\n";
                     $Sdai++;
                     next;
                 }

# "src/SdaiCONFIG_CONTROL_DESIGN.cc", line 7260: Warning: String literal converted to char* in formal argument 1 in call to SDAI::Select::Error(char*).

#  Filter responses over many lines 
#
#  SUN c++ warnings
#
#                print "$line";
                 if ( $line =~ m#In file included from\s+(.*):(\d+)([\,\:])(\s*)$# ) {
                     if ( $verbose ) {print "Start Multiline $line";}
                     chomp($line);
                     $continue=($3 =~m/\,/);
                     @loglines=($line);
                     unless ($continue) {print "Next line should have the message\n"};
                     while ($line=<GMAKE>) {
                         if ( $line =~ m#^\s+from\s+(.*):(\d+)([\,\:])(\s*)$# ) {
                             if ( $verbose ) {print "Continue Multiline $line";}
                             $continue=($3 =~m/\,/);
                         #   unless ($continue) {print "Next line should have the message\n"}
                             push(@loglines,$line);
                             next;
                         }
                         if ( $line =~ m#(.*):(\d+):\s+(\w+):(.*)$# ) {
                         #   print "The message $line ";
                             if ( $verbose ) {
                                 print "File $1\nLine $2\nSeverity$3\nWhinge $4\n";
                             }
                             print "WHINGE:  $4\n";
                             last;
                         }
                         print "NoMATCH terminating Multiline filter $line";
                         last;
                     }
                     next;
                 }
#
#  g3tog4 Fortran compilation (HP is verbose, linux runs ranlib)
#
#                print "$line";
                 if ( $line =~ m#libc stage done in# ) {
                     if ( $verbose ) {print "Start Fortran g3tog4: $line";}
                     chomp($line);
                     @f77lines=($line);
                     $f77error=0;
                     while ($line=<GMAKE>) {
                         chomp($line);
                         if ( $line =~ m#^\s*f77 # ) {
                             push(@f77lines,$line);
                             next;
                         }
                         if ( $line =~ m#^ar: Warning: creating# ) {
                             push(@f77lines,$line);
                             next;
                         }
                         if ( $line =~ m#^ar: creating# ) {
                             push(@f77lines,$line);
                             next;
                         }
                         if ( $line =~ m#^(/afs.*)\:$# ) {
                             $f77file=$1;
                             push(@f77lines,$line);
                             if ( $verbose ) {print "File $f77file for fort77 compilation $line";}
                             next;
                         }
                         if ( $line =~ m#^(/tmp/.*)\:$# ) {
                             $f77file=$1;
                             push(@f77lines,$line);
                             if ( $verbose ) {print "File $f77file for fort77 compilation $line";}
                             next;
                         }
                      #  if ( $line =~ m#^   (\w{4:12})\:# ) {
                      #  if ( $line =~ m#^   (\w{4:12})\:# ) {
                         if ( $line =~ m#^\s+(\w+)\:# ) {
                             $f77subprogram=$1;
                             if ( $verbose ) {print "Subprogram $f77subprogram";}
                             push(@f77lines,$line);
                             next;
                         }
                         if ( $line =~ m#^Creating/replacing object files in libG3toG4F.a$# ) {
                             push(@f77lines,$line);
                             next;
                         }
                         if ( $line =~ m#^Running ranlib on libG3toG4F.a$# ) {
                             push(@f77lines,$line);
                             next;
                         }
                         if ( $line =~ m#^libF stage done in# ) {
                             push(@f77lines,$line);
                             next;
                         }
                         if ( $line =~ m#^lib stage done in# ) {
                             push(@f77lines,$line);
                             if ( $f77error ) {
                                 print "G3toG4F $f77error problem lines in log\n";
                                 foreach $line (@f77lines) { print "  $line\n"; }
                             } else {
                                 print "G3toG4F completed with no errors logged.\n";
                             }
                             last;
                         }
                         print "Problem compiling G3toG4F: $line\n";
                         $f77error++;
                     }
                     next;
                 }
                 if ( $line =~ m#^Preprocessing.*\.ddl\s+\.\.\.# ) {
                    $line = &AnalyseSunddl($line,$filehandlemaybe);
                 }
                 print "$line\n";
#  Perhaps I can structure this code rather than just adding blocks.
#  need to call with the current line and maybe a filehandle to read
#  probably need to return the last line read to detect the end of my 
#  regime plus some "results"
             }
             close(GMAKE);
             print "ILD--------- $ild lines not printed\n";
         }
      }
   }
}
sub AnalyseSunddl {
    print "Text from Preprocessing  dll skipped\n";
    while ( $line = <GMAKE> ) { 
        if ( $line =~ /^gmake/ || $line =~ /Making dependency for file/ ) {return $line};
    #   print "ADDL: $line" 
    }
    print " AnalyseSundll read to EndOfFile\n";
}
                      

