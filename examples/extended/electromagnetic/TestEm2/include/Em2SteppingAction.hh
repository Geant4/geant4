// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em2SteppingAction.hh,v 1.2 1999-12-15 14:49:00 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em2SteppingAction_h
#define Em2SteppingAction_h 1

#include "G4UserSteppingAction.hh"

class Em2DetectorConstruction;
class Em2RunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em2SteppingAction : public G4UserSteppingAction
{
  public:
    Em2SteppingAction(Em2DetectorConstruction*,Em2RunAction*);
   ~Em2SteppingAction();

    void UserSteppingAction(const G4Step*);
  
  private:
  Em2DetectorConstruction* Em2Det;
  Em2RunAction*            Em2Run;  
};

#endif
