
#ifndef ExN04TrackerSD_h
#define ExN04TrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "ExN04TrackerHit.hh"
class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class ExN04TrackerSD : public G4VSensitiveDetector
{

  public:
      ExN04TrackerSD(G4String name);
      ~ExN04TrackerSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      ExN04TrackerHitsCollection *trackerCollection;

  public:
};




#endif

