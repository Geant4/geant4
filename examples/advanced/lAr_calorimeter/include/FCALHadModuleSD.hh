// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: FCALHadModuleSD.hh,v 1.3 2002-10-02 19:40:09 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef FCALHadModuleSD_h
#define FCALHadModuleSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class FCALTestbeamSetup;
class FCALHadModule;
class G4HCofThisEvent;
class G4Step;
#include "FCALCalorHit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class FCALHadModuleSD : public G4VSensitiveDetector
{
  public:
  
  //      FCALHadModuleSD(G4String, FCALTestbeamSetup* );
      FCALHadModuleSD(G4String); 
     ~FCALHadModuleSD();

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
      FCALHadModule* HadModule;
      G4int InitF2;
public:

  G4double EvisF2Tile[2330];

};

#endif

