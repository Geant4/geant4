//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef mySiSD_h
#define mySiSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class myDetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "mySensorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class mySiSD : public G4VSensitiveDetector
{
  public:
  
      mySiSD(G4String, myDetectorConstruction* );
     ~mySiSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
  
      mySensorHitsCollection*  SiCollection;   
      myDetectorConstruction* Si;
      G4int*                   HitSiID;
};

#endif





