#ifndef GammaRayTelPayloadSD_h
#define GammaRayTelPayloadSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class GammaRayTelDetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "GammaRayTelPayloadHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class GammaRayTelPayloadSD : public G4VSensitiveDetector
{
  public:
  
      GammaRayTelPayloadSD(G4String, GammaRayTelDetectorConstruction* );
     ~GammaRayTelPayloadSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
  
      GammaRayTelPayloadHitsCollection*  PayloadCollection;      
      GammaRayTelDetectorConstruction* Detector;
      G4int*                   HitID;
};

#endif




