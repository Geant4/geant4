#!/usr/local/bin/perl 
#

#print "TEST0\n";

$TestDir = "/afs/cern.ch/user/s/stesting/stt/dev2/Linux-g++/debug/stt/Linux-g++";

$flag_comp = 0;
$flag_warn = 0;

open(GMAKE,"$TestDir/gmake.log");

#print "TEST1\n";

while ($line=<GMAKE>) {

#    print "TEST2\n";
#    print $line;
    chomp($line);
    if ( $line =~ m#Compiling (\w+)\.# && $flag_comp eq 0) {
	 print "$line\n";
	 $flag_comp = 1;
	 $flag_warn = 0;
	 open(TMP,">>./buffer1.txt");
	 print TMP "$line\n";
	 close TMP;
	 next;
     }
    if ( $line =~ m#Compiling (\w+)\.# && $flag_comp eq 1 && $flag_warn eq 1) {
	 $flag_warn = 0;
	 open(TMP,">>./buffer2.txt");
	 close TMP;
	 system("cat ./buffer1.txt >> ./buffer2.txt");
	 open(TMP,">./buffer1.txt");
	 print TMP "$line\n";
	 close TMP;
	 next;
     }
    if ( $line !~ m#Compiling (\w+)\.# && $line !~ m#Creating/replacing object files in .*\/(\w+\.a)\s*# && $flag_comp eq 1) {
	 $flag_warn = 1;
	 $flag_lib  = 1;
	 open(TMP,">>./buffer1.txt");
	 print TMP "$line\n";
	 close TMP;
	 next;
     }
    if ( $line =~ m#Compiling (\w+)\.# && $flag_comp eq 1 && $flag_warn eq 0) {
	 open(TMP,">./buffer1.txt");
	 print TMP "$line\n";
	 close TMP;
	 next;
     }
    if ( $line =~ m#Creating/replacing object files in .*\/(\w+\.a)\s*# && $flag_lib eq 1) {
	 $flag_comp = 0;
	 $flag_warn = 0;
	 $flag_lib  = 0;
	 open(TMP,">>./buffer2.txt");
	 close TMP;
	 system("cat ./buffer1.txt >> ./buffer2.txt");
	 system("rm ./buffer1.txt");
	 system("mv ./buffer2.txt ./$1.txt");
	 next;
     }
    if ( $line =~ m#Creating/replacing object files in .*\/(\w+\.a)\s*# && $flag_lib eq 0) {
	 $flag_comp = 0;
	 $flag_warn = 0;
	 $flag_lib  = 0;
	 system("rm ./buffer1.txt");
	 next;
     }
    if ( $line !~ m#Compiling (\w+)\.# && $line !~ m#Creating/replacing object files in .*\/(\w+\.a)\s*# && $flag_comp eq 0) {
	 next;
     }
    
    
}
