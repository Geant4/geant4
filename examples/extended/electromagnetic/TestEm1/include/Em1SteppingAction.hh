// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1SteppingAction.hh,v 1.2 1999-12-15 14:48:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em1SteppingAction_h
#define Em1SteppingAction_h 1

#include "G4UserSteppingAction.hh"

class Em1RunAction;
class Em1EventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em1SteppingAction : public G4UserSteppingAction
{
  public:
    Em1SteppingAction(Em1RunAction*,Em1EventAction*);
   ~Em1SteppingAction();

    void UserSteppingAction(const G4Step*);
    
  private:
    Em1RunAction*   runAction;
    Em1EventAction* eventAction;    
};

#endif
