// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ZIIICalorimeterSD.hh,v 1.1 2001-06-26 11:24:07 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef ZIIICalorimeterSD_h
#define ZIIICalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"
#include "ZIIICalorHit.hh"

class ZIIIDetectorConstruction;
class G4Step;
class G4HCofThisEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class ZIIICalorimeterSD : public G4VSensitiveDetector
{
  public:
  
      ZIIICalorimeterSD(G4String);
     ~ZIIICalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*, G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
        void clear();
        void DrawAll();
        void PrintAll();

  private:
  
      ZIIICalorHitsCollection*  calorimeterCollection;      
};

#endif

