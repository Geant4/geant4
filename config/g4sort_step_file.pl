#!/usr/local/bin/perl

sub usage {
  print STDERR "This script takes a STEP file (from STDIN)\n";
  print STDERR "and writes it to STDOUT with the data section\n";
  print STDERR "sorted by ID.\n";
  exit 1;
}

sub error {
  print STDERR $_[0]."\n";
  exit 2;
}

if ($#ARGV >= 0) {
  &usage();
}

##########################################################################

$entry= "";
while (<STDIN>) {
  while ($_ ne "") {

    if (/[;']/) {            #']){

      $entry.= $`.$&;
      $_= $';
      if ($& eq ";") {

	# end of entry
	$lof= $_;
	&processEntry($entry);
	$entry= "";
	$_= $lof;

      } else {

	# start of string
	/[']/ || &error("String does not end within the line!"); #'];
	$entry.= $`.$&;
	$_= $';

      }
    } else {

      $entry.= $_;
      $_= "";

    }
  }
}

/^[ \t\n]*$/m || &error("File does not end with ';' !");
print "$_\n";
exit 0;

##########################################################################

sub processEntry {
  local($e)= @_;

  if (!$inData) {

    print $e;
    if ($e =~ /^(\s|\n)*DATA(\s|\n)*;$/m) {
      $inData= 1;
    }

  } else {

    if ($e =~ /^(\s|\n)*\#([0-9]+)\=/m) {

      $id= $2;
      $entries[$id]= $e;

    } else {

      $e =~ /^(\s|\n)*ENDSEC(\s|\n)*;$/m ||
	&error("Data section does not end with ENDSEC;");
      $inData= 0;
      for ($id=0; $id<=$#entries; $id++) {
	$en= $entries[$id];
	if ($en ne "") {
	  print $en;
	}
      }
      print $e;

    }
  }
}

##########################################################################
##########################################################################
##########################################################################


