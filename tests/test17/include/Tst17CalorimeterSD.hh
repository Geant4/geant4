// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst17CalorimeterSD.hh,v 1.1 1999-11-30 18:01:50 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst17CalorimeterSD_h
#define Tst17CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class Tst17DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "Tst17CalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst17CalorimeterSD : public G4VSensitiveDetector
{
  public:
  
      Tst17CalorimeterSD(G4String, Tst17DetectorConstruction* );
     ~Tst17CalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void PrintAll();

  private:
  
      Tst17CalorHitsCollection*  CalCollection;      
      Tst17DetectorConstruction* Detector;
      G4int*                   HitID;
};

#endif

