// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN03CalorimeterSD.hh,v 1.2 1999-12-15 14:49:23 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef ExN03CalorimeterSD_h
#define ExN03CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class ExN03DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "ExN03CalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class ExN03CalorimeterSD : public G4VSensitiveDetector
{
  public:
  
      ExN03CalorimeterSD(G4String, ExN03DetectorConstruction* );
     ~ExN03CalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
  
      ExN03CalorHitsCollection*  CalCollection;      
      ExN03DetectorConstruction* Detector;
      G4int*                   HitID;
};

#endif

