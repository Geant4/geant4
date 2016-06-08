#!/usr/bin/perl 

$filename = $ARGV[0];
open(INFILE, $filename) || die "\nSTEP file ".$filename." not found.";
print "Opening file ".$filename;
@filebuf = <INFILE>;
close(INFILE);

$g4 = "G4";
$outputfile = $g4 . $filename;
open(OUTFILE, ">$outputfile");
print "\nModifying STEP file for NIST toolkit...";

foreach  (@filebuf)
{
if(/.*VERTEX\_POINT.*/)
{
    ($a,$b) = split (",",$&);
    s//$a\,\'\'\,$b\n/;
}
elsif(/.*EDGE\_LOOP.*/)
{
    s/EDGE\_LOOP\(\'(\w+)\'\,/EDGE\_LOOP\(\'$1\'\,\'\'\,/;
}
elsif(/.*POLY\_LOOP.*/)
{
    s/POLY\_LOOP\(\'(\w+)\'\,/POLY\_LOOP\(\'$1\'\,\'\'\,/;
}
elsif(/.*ADVANCED\_FACE.*/)
{
    ($a,$b,$c) = split('\)',$&);
    s//$a\)\,\'\'$b\)$c\n/;
}
elsif(/.*FACE\_SURFACE.*/)
{
    ($a,$b,$c) = split('\)',$&);
    s//$a\)\,\'\'$b\)$c\n/;
}
elsif(/.*EDGE\_CURVE.*/)
{
    ($a,$b,$c,$d,$e) = split(",",$&);
    s//$a\,$b\,$c\,\'\'\,$d\,$e\n/;
}

print OUTFILE;
}

close(OUTFILE);
print "\nOutput written to ".$outputfile."\n";


