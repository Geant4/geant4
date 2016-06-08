#!/bin/csh -f
set file = $1
set tmpfile = fadsmv.tmp.$$
#if ( ! -f $file ) then
#  echo "File $file does not exist.  Stop."
#  exit 1
#endif

set filename = `basename $file`
set filepath = `dirname  $file`

echo $filename | sed \
  -e 's/DigitRootIO/G4DigitRootIO/g' \
  -e 's/EventRootIO/G4EventRootIO/g' \
  -e 's/HepMCRootIO/G4HepMCRootIO/g' \
  -e 's/HitRootIO/G4HitRootIO/g' \
  -e 's/MCTruthRootIO/G4MCTruthRootIO/g' \
  -e 's/RootFADSevent/G4RootEvent/g' \
  -e 's/RootIOManager/G4RootIOManager/g' \
  -e 's/RootTransManager/G4RootTransManager/g' \
  -e 's/fadsROOTLinkDef/G4ROOTLinkDef/g' \
-e 's/G4G4/G4/g' \
-e 's/G4VG4/G4V/g' \
  > $tmpfile

set filename2 = `cat $tmpfile`

if ( $filename != $filename2 ) then
  echo "mv -f $filepath/$filename $filepath/$filename2"
  mv -f $filepath/$filename $filepath/$filename2
endif

rm -f $tmpfile

