// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
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

#ifndef Em6SteppingAction_h
#define Em6SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"

class Em6DetectorConstruction;
class Em6RunAction;
class Em6EventAction;
class Em6SteppingMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em6SteppingAction : public G4UserSteppingAction
{
public: // Without description

    Em6SteppingAction(Em6DetectorConstruction*, Em6EventAction*,
                      Em6RunAction* );
   ~Em6SteppingAction();

    void UserSteppingAction(const G4Step*);

  private:
    Em6DetectorConstruction* detector;
    Em6EventAction*          eventaction;
    Em6RunAction*            runaction;
    Em6SteppingMessenger*    steppingMessenger;

    G4int IDnow,IDold;
    G4int evnoold ;

};

#endif
