// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst14CalorimeterSD.hh,v 1.3 1999-06-14 23:25:25 aforti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Tst14CalorimeterSD_h
#define Tst14CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class Tst14DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "Tst14CalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Tst14CalorimeterSD : public G4VSensitiveDetector
{
  public:
  
      Tst14CalorimeterSD(G4String, Tst14DetectorConstruction* );
     ~Tst14CalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void PrintAll();

  private:
  
      Tst14CalorHitsCollection*  CalCollection;      
      Tst14DetectorConstruction* Detector;
      G4int*                   HitID;
};

#endif

