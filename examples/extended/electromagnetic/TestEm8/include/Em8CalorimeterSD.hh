// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em8CalorimeterSD.hh,v 1.2 2000-06-27 13:29:46 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em8CalorimeterSD_h
#define Em8CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class Em8DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "Em8CalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em8CalorimeterSD : public G4VSensitiveDetector
{
  public:
  
      Em8CalorimeterSD(G4String, Em8DetectorConstruction* );
     ~Em8CalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void PrintAll();

  private:
  
      Em8CalorHitsCollection*  CalCollection;      
      Em8DetectorConstruction* Detector;
      G4int*                   HitID;
};

#endif

