// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: AnaEx01CalorimeterSD.hh,v 1.1.1.1 2000-09-14 11:37:21 barrand Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef AnaEx01CalorimeterSD_h
#define AnaEx01CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class AnaEx01DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "AnaEx01CalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class AnaEx01CalorimeterSD : public G4VSensitiveDetector
{
  public:
  
      AnaEx01CalorimeterSD(G4String, AnaEx01DetectorConstruction* );
     ~AnaEx01CalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
  
      AnaEx01CalorHitsCollection*  CalCollection;      
      AnaEx01DetectorConstruction* Detector;
      G4int*                   HitID;
};

#endif

