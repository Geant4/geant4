#ifndef DetectorSliceSensitiveEmCalo_h
#define DetectorSliceSensitiveEmCalo_h 1

#include "G4VSensitiveDetector.hh"
#include "DetectorSliceCalorimeterHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;


class DetectorSliceSensitiveEmCalo : public G4VSensitiveDetector {

public:

  DetectorSliceSensitiveEmCalo( const G4String name );
  ~DetectorSliceSensitiveEmCalo();

  void Initialize( G4HCofThisEvent* HCE );
  G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  void EndOfEvent( G4HCofThisEvent* HCE );

  void clear();
  void DrawAll();
  void PrintAll();

private:

  DetectorSliceCalorimeterHitsCollection* calCollection;
};


#endif
