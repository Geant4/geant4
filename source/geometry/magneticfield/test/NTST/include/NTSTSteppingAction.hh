// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: NTSTSteppingAction.hh,v 1.1 2003-11-07 21:30:28 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef NTSTSteppingAction_h
#define NTSTSteppingAction_h 1

#include "G4UserSteppingAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class NTSTSteppingAction : public G4UserSteppingAction
{
  public:
    NTSTSteppingAction();
   ~NTSTSteppingAction();

    void UserSteppingAction(const G4Step* aStep);
};

#endif
