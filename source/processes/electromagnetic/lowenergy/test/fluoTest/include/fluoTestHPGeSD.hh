//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef fluoTestHpGe_h
#define fluoTestHPGeSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class fluoTestDetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "fluoTestSensorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class fluoTestHPGeSD : public G4VSensitiveDetector
{
  public:
  
      fluoTestHPGeSD(G4String, fluoTestDetectorConstruction* );
     ~fluoTestHPGeSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
  
      fluoTestSensorHitsCollection*  HPGeCollection;   
      fluoTestDetectorConstruction* HPGe;
      G4int*                   HitHPGeID;
};

#endif





