#!/usr/local/bin/perl
#
#  Created:  Stefano Magni,  August 1998
# 
#   Modified: J.Apostolakis  September 11, 1998
#

# replaceIDs file

while (@ARGV){
   $file= shift @ARGV;

   open(FILE, $file) || die "Cannot open $file!\n";
   $tmp= "/tmp/lof";
   open(TMP, ">$tmp") || die "Cannot open temp file!\n";

   while (<FILE>) {

      $_ =~ s/\bTst10/Sc01/g ;
      print TMP $_;
   }
   close(TMP);
   close(FILE);
   system("mv -f $file .bak.$file");
   system("mv -f $tmp $file");
}
