//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoHpGe_h
#define XrayFluoHPGeSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"
#include "XrayFluoSensorHit.hh"

class XrayFluoDetectorConstruction;
class G4HCofThisEvent;
class G4Step;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoHPGeSD : public G4VSensitiveDetector
{
  public:
  
      XrayFluoHPGeSD(G4String, XrayFluoDetectorConstruction* );
     ~XrayFluoHPGeSD();

      void Initialize(G4HCofThisEvent*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

protected:

 G4bool ProcessHits(G4Step*,G4TouchableHistory*);
private:
  
      XrayFluoSensorHitsCollection*  HPGeCollection;   
      XrayFluoDetectorConstruction* Detector;
      G4int*                   HitHPGeID;
};

#endif





