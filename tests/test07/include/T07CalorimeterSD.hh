// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T07CalorimeterSD.hh,v 1.1 1999-01-08 16:35:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef T07CalorimeterSD_h
#define T07CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class T07DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "T07CalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class T07CalorimeterSD : public G4VSensitiveDetector
{
  public:
  
      T07CalorimeterSD(G4String, T07DetectorConstruction* );
     ~T07CalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
  
      T07CalorHitsCollection*  CalCollection;      
      T07DetectorConstruction* Detector;
      G4int*                   HitID;
};

#endif

