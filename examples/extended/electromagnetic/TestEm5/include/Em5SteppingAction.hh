// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em5SteppingAction.hh,v 1.2 1999-12-15 14:49:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em5SteppingAction_h
#define Em5SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"

class Em5DetectorConstruction;
class Em5RunAction;
class Em5EventAction;
class Em5SteppingMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em5SteppingAction : public G4UserSteppingAction
{
  public:
    Em5SteppingAction(Em5DetectorConstruction*, Em5EventAction*,
                      Em5RunAction* );
   ~Em5SteppingAction();

    void UserSteppingAction(const G4Step*);

  private:
    Em5DetectorConstruction* detector;
    Em5EventAction*          eventaction;
    Em5RunAction*            runaction;
    Em5SteppingMessenger*    steppingMessenger;

    G4int IDnow,IDold;
    G4int evnoold ;

};

#endif
