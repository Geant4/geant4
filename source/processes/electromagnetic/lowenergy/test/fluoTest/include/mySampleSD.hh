//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef mySampleSD_h
#define mySsmpleSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class myDetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "mySensorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class mySampleSD : public G4VSensitiveDetector
{
  public:
  
      mySampleSD(G4String, myDetectorConstruction* );
     ~mySampleSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();
 
  private:
  
      mySensorHitsCollection*  SamCollection;   
      myDetectorConstruction* Sample;
      G4int*                   HitSID;
};

#endif
