// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Class Description:
// The sensitive detector is defined
// Class Description - end
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Test17CalorimeterSD_h
#define Test17CalorimeterSD_h 1

#include "globals.hh"
#include "G4VSensitiveDetector.hh"

class Test17DetectorConstruction;
class G4HCofThisEvent;
class G4Step;

#include "Test17CalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test17CalorimeterSD : public G4VSensitiveDetector
{
public: // Without description
  
      Test17CalorimeterSD(G4String, Test17DetectorConstruction* );
     ~Test17CalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void PrintAll();

  private:
  
      Test17DetectorConstruction* detector;
};

#endif





