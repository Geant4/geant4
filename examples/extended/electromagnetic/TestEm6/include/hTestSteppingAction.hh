// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Class Description:
// Actions on each step are defined
// Class Description - end
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef hTestSteppingAction_h
#define hTestSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"

class hTestDetectorConstruction;
class hTestRunAction;
class hTestEventAction;
class hTestSteppingMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class hTestSteppingAction : public G4UserSteppingAction
{
public: // Without description

    hTestSteppingAction(hTestDetectorConstruction*, hTestEventAction*,
                      hTestRunAction* );
   ~hTestSteppingAction();

    void UserSteppingAction(const G4Step*);

  private:
    hTestDetectorConstruction* detector;
    hTestEventAction*          eventaction;
    hTestRunAction*            runaction;
    hTestSteppingMessenger*    steppingMessenger;

    G4int IDnow,IDold;
    G4int evnoold ;

};

#endif
