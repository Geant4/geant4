//
#ifndef BEAMTESTSILICONMONITOR_HH
#define BEAMTESTSILICONMONITOR_HH

#include "G4VSensitiveDetector.hh"
#include "BeamTestSiliconMonitorHit.hh"

class G4HCofThisEvent;
class G4TouchableHistory;
class G4Step;

class BeamTestSiliconMonitor : public G4VSensitiveDetector {

public:

  // Constructor
  BeamTestSiliconMonitor(const G4String& name);

  // Destructor
  virtual ~BeamTestSiliconMonitor();
  
  // Methods
  virtual void Initialize(G4HCofThisEvent* hitsCollectionOfThisEvent);
  
  virtual G4bool ProcessHits(G4Step* aStep,G4TouchableHistory* history);

  virtual void EndOfEvent(G4HCofThisEvent* hitsCollectionOfThisEvent);
  
private:
  
  // Data members
  BeamTestSiliconMonitorHitsCollection* fHitsCollection;
  G4int fHitsCollectionID;

};



#endif

