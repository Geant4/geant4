//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoSensorSD_h
#define  XrayFluoSensorSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class  XrayFluoDetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "XrayFluoSensorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class  XrayFluoSensorSD : public G4VSensitiveDetector
{
public:
  
  XrayFluoSensorSD(G4String,XrayFluoDetectorConstruction* );
  ~ XrayFluoSensorSD();
  
  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*,G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);
  void clear();
  void DrawAll();
  void PrintAll();
  
private:
  
  XrayFluoSensorHitsCollection*  SenCollection;   
  XrayFluoDetectorConstruction* Detector;
  G4int*                   HitID;
};

#endif





