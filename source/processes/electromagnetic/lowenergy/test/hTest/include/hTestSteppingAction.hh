#ifndef hTestSteppingAction_h
#define hTestSteppingAction_h 1

// -------------------------------------------------------------
//
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------
//      GEANT4 hTest
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- hTestSteppingAction -------------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of hTest 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4UserSteppingAction.hh"
#include "hTestDetectorConstruction.hh"
#include "hTestRunAction.hh"
#include "hTestEventAction.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestSteppingAction : public G4UserSteppingAction
{
public: // Without description

    hTestSteppingAction(hTestRunAction*, hTestEventAction*, 
                        hTestDetectorConstruction*);
   ~hTestSteppingAction();

    void UserSteppingAction(const G4Step*);

  private:
    hTestDetectorConstruction* detector;
    hTestEventAction*          eventaction;
    hTestRunAction*            runaction;

    G4int trIDnow;
    G4int trIDold;
    G4int evnOld;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

