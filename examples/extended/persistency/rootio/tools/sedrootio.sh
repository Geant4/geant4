#!/bin/csh -f
set file = $1
set tmpfile = $1.tmp
if ( ! -f $file ) then
  echo "File $file does not exist.  Stop."
  exit 1
endif

cat $file | sed \
-e 's/g4rootio/fadsrootio/g' \
-e 's/G4DCIOcatalog/DCIOcatalog/g' \
-e 's/G4DCIOentryT/DCIOentryT/g' \
-e 's/G4Pevent/FADSevent/g' \
-e 's/G4FileUtilities/FileUtilities/g' \
-e 's/G4HCIOcatalog/HCIOcatalog/g' \
-e 's/G4HCIOentryT/HCIOentryT/g' \
-e 's/G4PersistencyCenter/PersistencyCenter/g' \
-e 's/G4PersistencyCenterMessenger/PersistencyCenterMessenger/g' \
-e 's/G4PersistencyManager/PersistencyManager/g' \
-e 's/G4PersistencyManagerT/PersistencyManagerT/g' \
-e 's/G4VDCIOentry/VDCIOentry/g' \
-e 's/G4VHCIOentry/VHCIOentry/g' \
-e 's/G4VHepMCIO/VHepMCIO/g' \
-e 's/G4VMCTruthIO/VMCTruthIO/g' \
-e 's/G4VPDigitIO/VPDigitIO/g' \
-e 's/G4VPDigitsCollectionIO/VPDigitsCollectionIO/g' \
-e 's/G4VPEventIO/VPEventIO/g' \
-e 's/G4VPHitIO/VPHitIO/g' \
-e 's/G4VPHitsCollectionIO/VPHitsCollectionIO/g' \
-e 's/G4VTransactionManager/VTransactionManager/g' \
-e 's/G4MCTEvent/MCTEvent/g' \
-e 's/G4MCTEventAction/MCTEventAction/g' \
-e 's/G4MCTGenEvent/MCTGenEvent/g' \
-e 's/G4MCTGenParticle/MCTGenParticle/g' \
-e 's/G4MCTManager/MCTManager/g' \
-e 's/G4MCTMessenger/MCTMessenger/g' \
-e 's/G4MCTSimEvent/MCTSimEvent/g' \
-e 's/G4MCTSimParticle/MCTSimParticle/g' \
-e 's/G4MCTSimVertex/MCTSimVertex/g' \
-e 's/G4MCTSteppingAction/MCTSteppingAction/g' \
-e 's/G4MCTVStoreSelection/MCTVStoreSelection/g' \
-e 's/G4MCTVUserFilter/MCTVUserFilter/g' \
-e 's/G4HitRootIO/HitRootIO/g' \
-e 's/G4RootIOManager/RootIOManager/g' \
-e 's/G4RootTransManager/RootTransManager/g' \
-e 's/GEANT4 tag/FADS\/Goofy tag/g' \
-e 's/Geant4/FADS/g' \
-e 's/G4cout/std::cout/g' \
-e 's/G4cerr/std::cerr/g' \
-e 's/G4endl/std::endl/g' \
-e 's/G4std/std/g' \
-e 's/$package_name       = "FADS";/$package_name       = "FADS\/Goofy";/g' \
> $tmpfile

diff $file $tmpfile

