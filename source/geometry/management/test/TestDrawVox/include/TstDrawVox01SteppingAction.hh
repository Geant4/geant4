// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TstDrawVox01SteppingAction.hh,v 1.2 1999-12-15 14:49:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef TstDrawVox01SteppingAction_h
#define TstDrawVox01SteppingAction_h 1

#include "math.h"
#include "G4UserSteppingAction.hh"
#include "Histo.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class TstDrawVox01SteppingAction : public G4UserSteppingAction
{
  public:
    TstDrawVox01SteppingAction();
    virtual ~TstDrawVox01SteppingAction();

    virtual void UserSteppingAction(const G4Step*);
private:
  odHisto Steplength;
  odHisto SteplengthProfile;
};

#endif
