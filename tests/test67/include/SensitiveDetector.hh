#ifndef Sensitive_h
#define Sensitive_h 1

#include "G4VSensitiveDetector.hh"
#include "Hits.hh"

class G4Step;
class G4HCofThisEvent;

class SensitiveDetector : public G4VSensitiveDetector
{
  public:
      SensitiveDetector(G4String);
     ~SensitiveDetector();

      void Initialize(G4HCofThisEvent*);
  
      G4bool ProcessHits(G4Step*, G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);

  private:
      HitsCollection *Collection;
  

};

#endif
