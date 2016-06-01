
#ifndef ExN04MuonSD_h
#define ExN04MuonSD_h 1

#include "G4VSensitiveDetector.hh"
#include "ExN04MuonHit.hh"
class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class ExN04MuonSD : public G4VSensitiveDetector
{

  public:
      ExN04MuonSD(G4String name);
      ~ExN04MuonSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      ExN04MuonHitsCollection * muonCollection;
      G4double positionResolution;

};




#endif

