#! /usr/bin/perl

while(<STDIN>) {
  chop;
  s/.hh//;
  print "  -e \'s/$_/G4$_/g\' \\\n";

}

