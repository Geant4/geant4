//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef fluoTestSiSD_h
#define fluoTestSiSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class fluoTestDetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "fluoTestSensorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class fluoTestSiSD : public G4VSensitiveDetector
{
  public:
  
      fluoTestSiSD(G4String, fluoTestDetectorConstruction* );
     ~fluoTestSiSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
  
      fluoTestSensorHitsCollection*  SiCollection;   
      fluoTestDetectorConstruction* Si;
      G4int*                   HitSiID;
};

#endif





