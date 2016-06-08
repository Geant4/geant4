#!/bin/csh -f
set file = $1
set tmpfile = $1.tmp.$$
if ( ! -f $file ) then
  echo "File $file does not exist.  Stop."
  exit 1
endif

cat $file | sed \
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
-e 's/FADS\/Goofy tag/GEANT4 tag/g' \
  > $tmpfile

diff $file $tmpfile

mv -f $file ${file}.org
mv -f $tmpfile $file

