//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FluoTestHpGe_h
#define FluoTestHPGeSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class FluoTestDetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "FluoTestSensorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FluoTestHPGeSD : public G4VSensitiveDetector
{
  public:
  
      FluoTestHPGeSD(G4String, FluoTestDetectorConstruction* );
     ~FluoTestHPGeSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
  
      FluoTestSensorHitsCollection*  HPGeCollection;   
      FluoTestDetectorConstruction* Detector;
      G4int*                   HitHPGeID;
};

#endif





