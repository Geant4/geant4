//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef XrayFluoSampleSD_h
#define XrayFluoSampleSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class XrayFluoDetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "XrayFluoSensorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class XrayFluoSampleSD : public G4VSensitiveDetector
{
public:
  
  XrayFluoSampleSD(G4String, XrayFluoDetectorConstruction* );
  ~XrayFluoSampleSD();
  
  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*,G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);
  void clear();
  void DrawAll();
  void PrintAll();
  
private:
  
  XrayFluoSensorHitsCollection*  SamCollection;   
  XrayFluoDetectorConstruction* Sample;
  G4int*                   HitSID;
};

#endif
