// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em10CalorimeterSD.hh,v 1.1 2000-07-14 15:51:14 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em10CalorimeterSD_h
#define Em10CalorimeterSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class Em10DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "Em10CalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em10CalorimeterSD : public G4VSensitiveDetector
{
  public:
  
      Em10CalorimeterSD(G4String, Em10DetectorConstruction* );
     ~Em10CalorimeterSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void PrintAll();

  private:
  
      Em10CalorHitsCollection*  CalCollection;      
      Em10DetectorConstruction* Detector;
      G4int*                   HitID;
};

#endif

