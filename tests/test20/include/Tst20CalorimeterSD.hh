// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20CalorimeterSD.hh,v 1.1 2001-05-25 12:50:05 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst20CalorimeterSD_h
#define Tst20CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class Tst20DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "Tst20CalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst20CalorimeterSD : public G4VSensitiveDetector
{
  public:
  
      Tst20CalorimeterSD(G4String, Tst20DetectorConstruction* );
     ~Tst20CalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void PrintAll();

  private:
  
      Tst20CalorHitsCollection*  CalCollection;      
      Tst20DetectorConstruction* Detector;
      G4int*                   HitID;
};

#endif

