
#ifndef ExN04CalorimeterSD_h
#define ExN04CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "ExN04CalorimeterHit.hh"
class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class ExN04CalorimeterSD : public G4VSensitiveDetector
{

  public:
      ExN04CalorimeterSD(G4String name);
      ~ExN04CalorimeterSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
      ExN04CalorimeterHitsCollection *CalCollection;
      int CellID[20][48];
      const int numberOfCellsInZ;
      const int numberOfCellsInPhi;
};




#endif

