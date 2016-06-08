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
    s/VERTEX\_POINT\s*\(\s*\'\s*\'\s*\,/VERTEX\_POINT\(\'\'\,\'\'\,/g ;
    s/VERTEX\_POINT\s*\(\s*\'NONE\'\s*\,/VERTEX\_POINT\(\'NONE\'\,\'\'\,/g ;
}
elsif(/.*EDGE\_LOOP.*/)
{
    s/EDGE\_LOOP\s*\(\s*\'\s*\'\s*\,/EDGE\_LOOP\(\'\'\,\'\'\,/g;
    s/EDGE\_LOOP\s*\(\s*\'NONE\'\s*\,/EDGE\_LOOP\(\'NONE\'\,\'\'\,/g;
}
elsif(/.*POLY\_LOOP.*/)
{
    s/POLY\_LOOP\s*\(\s*\'\s*\'\s*\,/POLY\_LOOP\(\'\'\,\'\'\,/g;
    s/POLY\_LOOP\s*\(\s*\'NONE\'\s*\,/POLY\_LOOP\(\'NONE\'\,\'\'\,/g;
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


