// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TstVASteppingAction.cc,v 1.2 2001-02-01 21:27:20 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "TstVASteppingAction.hh"
#include "G4SteppingManager.hh"
#include "math.h"
#include "g4std/fstream"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TstVASteppingAction::TstVASteppingAction() : Steplength(100,0.,100.),
SteplengthProfile(100,0.,2*M_PI)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TstVASteppingAction::~TstVASteppingAction()
{
  G4std::ofstream o("test01.stepLength.plt");
  Steplength.output(o);
  o.close();
  o.open("test01.stepLengthProfile.plt");
  SteplengthProfile.output(o);
  o.close();
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TstVASteppingAction::UserSteppingAction(const G4Step* aStep)
{
  Steplength.accumulate(aStep->GetStepLength());

  G4double phi = aStep->GetDeltaPosition().phi();
  if (phi < 0.) phi = phi + twopi;
  SteplengthProfile.accumulate(phi,aStep->GetStepLength());

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....















