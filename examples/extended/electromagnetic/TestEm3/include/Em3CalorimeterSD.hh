// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em3CalorimeterSD.hh,v 1.2 1999-12-15 14:49:01 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em3CalorimeterSD_h
#define Em3CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class Em3DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "Em3CalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em3CalorimeterSD : public G4VSensitiveDetector
{
  public:
  
      Em3CalorimeterSD(G4String, Em3DetectorConstruction* );
     ~Em3CalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
  
      Em3CalorHitsCollection*  CalCollection;      
      Em3DetectorConstruction* Detector;
      G4int*                   HitID;
};

#endif

