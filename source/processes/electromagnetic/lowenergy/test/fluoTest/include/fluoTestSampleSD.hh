//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef fluoTestSampleSD_h
#define fluoTestSsmpleSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class fluoTestDetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "fluoTestSensorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class fluoTestSampleSD : public G4VSensitiveDetector
{
  public:
  
      fluoTestSampleSD(G4String, fluoTestDetectorConstruction* );
     ~fluoTestSampleSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();
 
  private:
  
      fluoTestSensorHitsCollection*  SamCollection;   
      fluoTestDetectorConstruction* Sample;
      G4int*                   HitSID;
};

#endif
