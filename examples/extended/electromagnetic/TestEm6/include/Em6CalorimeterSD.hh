// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em6CalorimeterSD_h
#define Em6CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class Em6DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "Em6CalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em6CalorimeterSD : public G4VSensitiveDetector
{
  public:
  
      Em6CalorimeterSD(G4String, Em6DetectorConstruction* );
     ~Em6CalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void PrintAll();

  private:
  
      Em6CalorHitsCollection*  CalCollection;      
      Em6DetectorConstruction* Detector;
      G4int*                   HitID;
};

#endif

