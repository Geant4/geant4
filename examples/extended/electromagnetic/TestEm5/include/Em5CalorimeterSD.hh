// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em5CalorimeterSD.hh,v 1.1 1999-10-12 12:23:29 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em5CalorimeterSD_h
#define Em5CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class Em5DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "Em5CalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em5CalorimeterSD : public G4VSensitiveDetector
{
  public:
  
      Em5CalorimeterSD(G4String, Em5DetectorConstruction* );
     ~Em5CalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void PrintAll();

  private:
  
      Em5CalorHitsCollection*  CalCollection;      
      Em5DetectorConstruction* Detector;
      G4int*                   HitID;
};

#endif

