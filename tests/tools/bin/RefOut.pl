#!/usr/local/bin/perl
#
# Move 'reference' output to CVS working dir
# (you should commit it by hand [current!])

# SUN-CC 'reference' output in $G4WORKDIR will be moved to $G4INSTALL area 
# in 'prod' test area

$FromDir = "/afs/cern.ch/sw/geant4/stt/prod/SUN-CC/debug_NONISO/stt/SUN-CC";
$ToDir   ="/afs/cern.ch/sw/geant4/stt/prod/src/geant4/tests"; 


opendir(FROM,$FromDir);
@fromFiles = readdir(FROM);
closedir(FROM);


%map = ();
foreach $file (@fromFiles) {
#    print "$file\n";

    if ($file =~ m/^test[0-9]*(|\.hadron|\.EMtest)\.out/) {
#	print "$file\n";
	($dir,$others) = split(/\./,$file);
#	print "$dir\n";
	%map = ($file=>$dir,%map);
    }
}

while (($file,$dir) = each (%map)) {


    $real_file = $file;

    $link = readlink("$ToDir/$dir/$file");
    if ($link){
	print "link=$link\n";
	$real_file = $link;
    }

# 
# Update to HEAD (remove sticky tag)
    system("cvs update -A $ToDir/$dir/$real_file");
    system("cp $FromDir/$file $ToDir/$dir/$file");
    system("cvs ci -m 'New reference output.' $ToDir/$dir/$real_file");
#    system("ls $FromDir/$file");
#    system("ls $ToDir/$dir/$file");


#    opendir(TO,"$ToDir/$dir");
#    @toFiles = readdir(TO);
#    closedir(TO);
#    
#    foreach $tofile (@toFiles) {
#	print "$tofile\n";
#    }
    
#    print "$out\n";
    print "$file=$dir\n";
}

