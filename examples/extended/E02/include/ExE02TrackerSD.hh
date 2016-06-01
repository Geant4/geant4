
#ifndef ExE02TrackerSD_h
#define ExE02TrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "ExE02TrackerHit.hh"
#include "G4Step.hh"

class ExE02TrackerSD : public G4VSensitiveDetector
{

  public:
      ExE02TrackerSD(G4String name);
      ~ExE02TrackerSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      ExE02TrackerHitsCollection *EvenCollection;
      ExE02TrackerHitsCollection *OddCollection;
};


#endif

