#ifndef TARGETLAYERSD_HH
#define TARGETLAYERSD_HH

#include "TargetLayerHitsCollection.hh"
#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;


class TargetLayerSD : public G4VSensitiveDetector {

 public:
   TargetLayerSD(G4String);
   ~TargetLayerSD();

   void Initialize(G4HCofThisEvent*);
   void EndOfEvent(G4HCofThisEvent*);

 private:
   G4bool ProcessHits(G4Step*, G4TouchableHistory*);
   TargetLayerHitsCollection* hitsCollection;

   G4int hitsCollectionID;
}; 

#endif // TARGETLAYERSD_HH
