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

#ifndef Test17SteppingAction_h
#define Test17SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"

class Test17DetectorConstruction;
class Test17RunAction;
class Test17EventAction;
class Test17SteppingMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Test17SteppingAction : public G4UserSteppingAction
{
public: // Without description

    Test17SteppingAction(Test17DetectorConstruction*, Test17EventAction*,
                      Test17RunAction* );
   ~Test17SteppingAction();

    void UserSteppingAction(const G4Step*);

  private:
    Test17DetectorConstruction* detector;
    Test17EventAction*          eventaction;
    Test17RunAction*            runaction;
    Test17SteppingMessenger*    steppingMessenger;

    G4int IDnow,IDold;
    G4int evnoold ;

};

#endif
