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

#ifndef hTestCalorimeterSD_h
#define hTestCalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class hTestDetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "hTestCalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestCalorimeterSD : public G4VSensitiveDetector
{
public: // Without description
  
      hTestCalorimeterSD(G4String, hTestDetectorConstruction* );
     ~hTestCalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void PrintAll();

  private:
  
      hTestCalorHitsCollection*  CalCollection;      
      hTestDetectorConstruction* Detector;
      G4int HitID[500];
};

#endif

