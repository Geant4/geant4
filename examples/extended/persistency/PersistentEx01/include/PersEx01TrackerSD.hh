
#ifndef PersEx01TrackerSD_h
#define PersEx01TrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "PersEx01TrackerHit.hh"
#include "G4Step.hh"

class PersEx01TrackerSD : public G4VSensitiveDetector
{

  public:
      PersEx01TrackerSD(G4String name);
      ~PersEx01TrackerSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      PersEx01TrackerHitsCollection *EvenCollection;
      PersEx01TrackerHitsCollection *OddCollection;
};


#endif

