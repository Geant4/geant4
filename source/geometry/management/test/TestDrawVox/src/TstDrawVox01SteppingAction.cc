// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TstDrawVox01SteppingAction.cc,v 1.1 1999-07-28 17:57:27 graignac Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "TstDrawVox01SteppingAction.hh"
#include "G4SteppingManager.hh"
#include "math.h"
#include <fstream.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TstDrawVox01SteppingAction::TstDrawVox01SteppingAction() : Steplength(100,0.,100.),
SteplengthProfile(100,0.,2*M_PI)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TstDrawVox01SteppingAction::~TstDrawVox01SteppingAction()
{
  ofstream o("test01.stepLength.plt");
  Steplength.output(o);
  o.close();
  o.open("test01.stepLengthProfile.plt");
  SteplengthProfile.output(o);
  o.close();
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TstDrawVox01SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  Steplength.accumulate(aStep->GetStepLength());

  G4double phi = aStep->GetDeltaPosition().phi();
  if (phi < 0.) phi = phi + twopi;
  SteplengthProfile.accumulate(phi,aStep->GetStepLength());

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....















