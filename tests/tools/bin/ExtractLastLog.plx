#!/usr/local/bin/perl
#
# /afs/cern.ch/sw/geant4/stt/dev1/testtools/geant4/tests/tools/bin
# ExtractLastLog.plx
#
# 20-Nov-2000 stthosts.data table changed to configure for machines

$ActiveExamination="doit";
$ShowStt=1;

open(CONFIG,"OnTest") || die "Failed to open OnTest configuration file $! ";
($DevDir,$Tag)=split(' ',<CONFIG>);
close(CONFIG);
print "\nWorking in dir \"$DevDir\" with test-set tag \"$Tag\"\n";
$SttTag="stt." . "$Tag";
$TestTop="/afs/cern.ch/sw/geant4/stt/$DevDir";
$TestLogDir="/afs/cern.ch/sw/geant4/stt/$DevDir/testtools/geant4/tests/tools/bin";

&ReadConfigurationFiles($DevDir);          # @platforms  @tests and %HostName are global.
&SetInitialTable();

opendir(TL,"$TestLogDir") || die "Failed to opendir TestLog  $TestLogDir $!";
@testlogs=grep(m/^\w+\.(dev\d|prod)\.\w+.*\.log/,readdir(TL));
closedir(TL);
foreach $testlog (@testlogs) {
    next unless ((-M "$TestLogDir/$testlog") < 60 );
#   print "Test Run Log \"$testlog\"\n";
    $Machine="Machine";
    $Option="CompilerOpts";
    if ( $testlog =~ m/(\w+)\.\w+\.(.*)\.log/ ) {
        $Machine=$1;
        $Option=$2;
    }
    $lines=0;
    $title=0;
    $copy=0;
    undef($Step); undef($ResultsDir); 
    undef(%Start); undef(%Disabled); undef(%Finish); undef(%Missing);
    open(TLC,"$TestLogDir/$testlog") || die "Failed to open read $TestLogDir/$testlog $!";
    while ($line = <TLC> ) {
        $lines++;
        if ( $line =~ /^_____/ ) {$copy=0;
            # next; Problem with metacharacter in tag - use index on next line for now
            $newtestline=$line;
            $line = <TLC>;
            $copy = (index($line,$Tag) > 0);
        }
        if ( $copy ) {
            if ( 0 == $title ) { print "\n\nTag:  $Tag\nTest:  $testlog\n\n$newtestline";$title=2;}
            $lines++;
            chomp($line);
            if ( $ShowAll ) {print "$Machine $Option $lines    $line\n"};
            if ( $line =~ m#^G4WORKDIR.*\/(\S+\/\S+)$# ) {
                $platform=$1;
                $platform=~s#/#.#;
                &change_table("Hosts",$platform,$Machine);
            }
            if ( $line =~ m/^STT:(\w+)\s+(\w+)/ ) {
                if ( $ShowStt ) {print "$Machine $Option $lines    $line\n"}
                $Step=$1; $State=$2;
# Build for now means build test executables after library compilation
#        logged in gmake.log.
                if ( $line =~ m/^STT:BUILD Started/ ) {
                     $nextline = <TLC>;
                     &change_table("Compile",$platform,"S");
                     &change_table("Build",$platform,"S");
                     $lines++;
                     if ( $nextline =~ m#^(/afs/.*)\s+created# ) {
                         $ResultsDir=$1;
                         undef($syntaxerror);
                     } else {
                         print "Unexpected text in log file\n$nextline\ncreation of workdir expected\n";
                     }
                }
                if ( $line =~ m/^STT:BUILD Finished/ ) {
                     &change_table("Compile",$platform,"SF");
                     &change_table("Build",$platform,"SF");
                }
            }
            if ( $line =~ m/^In progress already/ ) {
                print "Do something neat here .... \n";
                $nextline = <TLC>;
                $lines++;
                print "This should be the ls output on the locking file\n$nextline";
            }
            if ( $line =~ m/^Starting (test\d+) in/ ) { $Start{$1}=$line }
            if ( $line =~ /Finished (test\d+) in/ ) { $Finish{$1}=$line }
            if ( $line =~ /Disabled (test\d+) in/ ) { $Disabled{$1}=$line }
            if ( $line =~ /Missing (test\d+) in/ ) { $Missing{$1}=$line }
            if ( $line =~ m/^Starting (test\d+\.\w+) in/ ) { $Start{$1}=$line }
            if ( $line =~ /Finished (test\d+\.\w+) in/ ) { $Finish{$1}=$line }
            if ( $line =~ /Disabled (test\d+\.\w+) in/ ) { $Disabled{$1}=$line }
            if ( $line =~ /Missing (test\d+\.\w+) in/ ) { $Missing{$1}=$line }
            if ( $line =~ m/syntax\s+error/ ) {
                $syntaxerror=$line;
                print "syntax error reported\n$line\n";
            }
        }
    }
    close(TLC);
    next unless(defined($Step));
    print "Reached $Step $State Results Directory is ...\n";
    print "$ResultsDir\n";
    print "\n";
    next if ( $Machine =~ /dxplus00/ );
    if ( $ActiveExamination && "$Step.$State" eq "BUILD.Started" ) {
        @MakingProgress=`rsh -l stesting $Machine "/afs/cern.ch/sw/geant4/stt/dev1/testtools/geant4/tests/tools/bin/FileAge.plx $ResultsDir/gmake.log"`;
        foreach $line (@MakingProgress) {
            chomp($line);
            if ( $line =~ m#(\d+)\s+(\d+)\s+/# ) {
                &change_table("Compile",$platform,"$1 $2");
                print "Gmake Age $line\n";
            }
        }
    }
    if (defined($syntaxerror)) {
        print "$syntaxerror\n";
        print "SE-- $ResultsDir\n";
    }
    undef($pl);undef($sl);
    foreach $testnum (keys (%Start)) {
        $TestState{$testnum}="S";
        if (defined($Finish{$testnum})) {$TestState{$testnum}.="F"}
        &change_table("$testnum",$platform,"$TestState{$testnum}");
        $num=$testnum;
        $num=~s/test//;
        $pl.=sprintf("%4d",$num);
        $sl.=sprintf("%4s",$TestState{$testnum});
    }
    foreach $testnum (keys (%Disabled)) {
        $TestState{$testnum}="D";
        &change_table("$testnum",$platform,"$TestState{$testnum}");
        $num=$testnum;
        $num=~s/test//;
        $pl.=sprintf("%4d",$num);
        $sl.=sprintf("%4s",$TestState{$testnum});
    }
    foreach $testnum (keys (%Missing)) {
        $TestState{$testnum}="M";
        &change_table("$testnum",$platform,"$TestState{$testnum}");
        $num=$testnum;
        $num=~s/test//;
        $pl.=sprintf("%4d",$num);
        $sl.=sprintf("%4s",$TestState{$testnum});
    }
    print "Test  $pl\nState $sl\n" if (defined($pl));
    
}
print "\n\n\n\n\n =================================================\n\n\n\n";
#  &print_all_entries_for_test();
&open_html($Tag,$DevDir);
&print_head($Tag,$DevDir);
&print_platform_headings_and_state("Hosts");
&print_state_for_test("Compile");
&print_state_for_test("Build");
&print_state_for_all_tests();
&print_tail($Tag,$DevDir);
&close_html($Tag,$DevDir);
exit();

sub ReadConfigurationFiles {
# Argument is $devdir (prod dev1 dev2)
    @tests=("Hosts","Compile","Build");
    open(TESTS,"stttests.data") || die "Failed to open (read) stttests.data $!";
    while ($line = <TESTS> ) {
        if ( $line =~ m/^(test.*)\s/ ) { push(@tests,$1) }
    }
    close(TESTS);

    open(HOSTS,"stthosts.data") || die "Failed to open (read) stthosts.data $!";
    while ($line = <HOSTS> ) {
        next if ($line =~ m/^#/);
        next if ($line =~ m/^%/);
        if ( $line =~ m/^$_[0]\s+(\S+)\s+(\S+)\s+(\S+)\s+/ ) { 
            print "$1 $2 $3 $4\n";
            push(@platforms,"$2.$3");
            $HostName{"$2.$3"}=$1;
        }
    }
    close(HOSTS);
}

foreach $platform (sort(@platforms)) { print "$platform\n";}

sub SetInitialTable{

    foreach $test (@tests) {
        foreach $platform (sort(@platforms)) { 
            &update_table($test,$platform,"--");
        }
    }

    $test="Hosts";
    foreach $platform (sort(@platforms)) { 
        &change_table($test,$platform,$HostName{$platform});
    }
}

# these are from the Panther book ?
# foreach $ymbj (@{$test_index{$test}}) {print "ymbj @{$ymbj} ymbj\n";}
# I don't understand them yet.

sub change_table {
    my($test,$platform,$state)=@_;
    $rlwork  = [$test,$platform,$state];
# I don't want to push a new value - I want to change the old one !
    foreach $rlwork (@{$test_index{$test}}) {
        if ( "$rlwork->[1]" eq "$platform" ) {
            $rlwork->[2]=$state;
        }
    }
# we the above works (from the test point of view) but I don't 
# think this is the "right" way to do it.
}

sub update_table {
    my($test,$platform,$state)=@_;
    $rlEntry = [$test,$platform,$state];
    push(@{$test_index{$test}}, $rlEntry);
    push(@{$platform_index{$platform}}, $rlEntry);
}

sub print_row_for_test {
    my($test)=@_;
    print ("TestStep: $test\n");
    foreach $rlEntry (@{$test_index{$test}}) {
        print("\t", $rlEntry->[1], " : ",$rlEntry->[2], "\n");
    }
}
sub print_state_for_test {
    my($test)=@_;my($l1);
    print HTML "<TR>\n";
    print HTML "<TD>$test</TD>\n";
    $l1 = sprintf("%-7s ",substr($test,0,7));
    foreach $rlEntry (@{$test_index{$test}}) {
        print HTML ("<TD>" ,$rlEntry->[2], "</TD>");
        $l1 .= sprintf("  %-6s",$rlEntry->[2]);
    }
    print HTML "</TR>\n";
    print "$l1\n";
}
sub print_platform_headings_and_state{
    my($test)=@_;
    my($architecture);my($compileroptions);my($hostname);my($pco);
    my($l1),my($l2),my($l3);
    print HTML "<TABLE>\n";
    print HTML "<TR><TD>TestStage</TD>\n";
    $l1="        ";
    $l2="  Test  ";
    $l3="        ";
    foreach $rlEntry (@{$test_index{$test}}) {
        ($architecture,$compileroptions)=split('\.',$rlEntry->[1]);
        $hostname=$rlEntry->[2];
        $pco=substr($compileroptions,0,1) . "-x"; # when (non)iso not specified
        if ($compileroptions =~ m/NON/ ) {$pco=~s/-x/-n/};
        if ($compileroptions =~ m/ISO/ ) {$pco=~s/-x/-i/};
        print HTML "<TD>$architecture $pco   $hostname</TD>";
        $l1 .= sprintf("%-7s ",substr($architecture,0,7));
        $l2 .= sprintf("  %-6s", $pco);
        $hostname=~s/plus/+/;
        $l3 .= sprintf("%-7s ",substr($hostname,0,7));
    }
    print HTML "</TR>\n";
    print "$l1\n$l2\n$l3\n";
}

sub print_all_entries_for_test {
    foreach $test ( sort (keys(%test_index)) ) {
        next unless( $test =~ m/^test/ );
        &print_row_for_test($test);
    }
}
sub print_state_for_all_tests {
    foreach $test ( sort (keys(%test_index)) ) {
        next unless( $test =~ m/^test/ );
        &print_state_for_test($test);
    }
}
sub print_head {
    my($Tag);my($DevDir);
    print "Progress On Test $_[0] in $_[1]\n";
    print HTML "<HTML><HEAD><TITLE>\n";
    print HTML "Geant4 STT Progress On Test $_[0] in $_[1]\n";
    print HTML "</TITLE></HEAD>\n<BODY>\n";
    print HTML "<H3>Progress On Test $_[0] in $_[1]</H3>\n";
    print HTML "<P>\n";
    print HTML "Check iso date formats\n";
    print HTML "</P>\n";
}
sub print_tail {
    print HTML "</BODY>\n<?HTML>\n";
}
sub open_html {
    $filename="$_[1].html";
    open(HTML,">$filename") || print "Failed to open write $filename $!\n";
}
sub close_html {
#   $filename="$_[1].html";
    close(HTML,">$filename"); # rename and move earlier copies next
}
