
#ifndef MyTrackerSD_h
#define MyTrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "MyTrackerHit.hh"
#include "G4Step.hh"

class MyTrackerSD : public G4VSensitiveDetector
{

  public:
      MyTrackerSD(G4String name);
      ~MyTrackerSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      MyTrackerHitsCollection *EvenCollection;
      MyTrackerHitsCollection *OddCollection;
};


#endif

