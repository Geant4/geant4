// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F03CalorimeterSD.hh,v 1.1 2001-06-08 11:55:39 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F03CalorimeterSD_h
#define F03CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class F03DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "F03CalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F03CalorimeterSD : public G4VSensitiveDetector
{
  public:
  
      F03CalorimeterSD(G4String, F03DetectorConstruction* );
     ~F03CalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void PrintAll();

  private:
  
      F03CalorHitsCollection*  CalCollection;      
      F03DetectorConstruction* Detector;
      G4int*                   HitID;
};

#endif

