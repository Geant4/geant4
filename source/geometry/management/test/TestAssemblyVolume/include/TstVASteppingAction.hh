// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TstVASteppingAction.hh,v 1.2 2001-02-01 21:25:40 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef TstVASteppingAction_h
#define TstVASteppingAction_h 1

#include "math.h"
#include "G4UserSteppingAction.hh"
#include "Histo.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class TstVASteppingAction : public G4UserSteppingAction
{
  public:
    TstVASteppingAction();
    virtual ~TstVASteppingAction();

    virtual void UserSteppingAction(const G4Step*);
private:
  odHisto Steplength;
  odHisto SteplengthProfile;
};

#endif
