// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em0SteppingAction.hh,v 1.2 1999-05-10 16:15:11 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em0SteppingAction_h
#define Em0SteppingAction_h 1

#include "G4UserSteppingAction.hh"

class Em0EventAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em0SteppingAction : public G4UserSteppingAction
{
  public:
    Em0SteppingAction(Em0EventAction*);
   ~Em0SteppingAction();

    void UserSteppingAction(const G4Step *);
    
  private:
    Em0EventAction* eventAction;    
};

#endif
