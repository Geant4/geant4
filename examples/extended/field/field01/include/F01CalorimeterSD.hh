// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: F01CalorimeterSD.hh,v 1.1 2001-03-27 16:21:24 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F01CalorimeterSD_h
#define F01CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class F01DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "F01CalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F01CalorimeterSD : public G4VSensitiveDetector
{
  public:
  
      F01CalorimeterSD(G4String, F01DetectorConstruction* );
     ~F01CalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void PrintAll();

  private:
  
      F01CalorHitsCollection*  CalCollection;      
      F01DetectorConstruction* Detector;
      G4int*                   HitID;
};

#endif

