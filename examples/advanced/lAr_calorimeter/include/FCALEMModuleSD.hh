// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FCALEMModuleSD.hh,v 1.3 2002-10-02 19:40:09 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FCALEMModuleSD_h
#define FCALEMModuleSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"
#include "G4Event.hh"

class FCALTestbeamSetup;
class FCALEMModule;
class G4HCofThisEvent;
class G4Step;
#include "FCALCalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FCALEMModuleSD : public G4VSensitiveDetector
{
  public:
  
  //      FCALEMModuleSD(G4String, FCALTestbeamSetup* );
      FCALEMModuleSD(G4String); 
     ~FCALEMModuleSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
  
  //  FCALCalorHitsCollection*  CalCollection;      
  //      FCALTestbeamSetup* Detector;
  //   G4int*                   HitID;
  FCALEMModule* EmModule;
  G4int Init_state;
  
  public:

  G4double EvisF1Tile[1131];

};

#endif

