// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TstVASteppingAction.cc,v 1.3 2001-02-07 17:31:02 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------

#include "TstVASteppingAction.hh"
#include "G4SteppingManager.hh"
#include "g4std/fstream"

TstVASteppingAction::TstVASteppingAction() : Steplength(100,0.,100.),
SteplengthProfile(100,0.,2*M_PI)
{ }

TstVASteppingAction::~TstVASteppingAction()
{
  G4std::ofstream o("test01.stepLength.plt");
  Steplength.output(o);
  o.close();
  o.open("test01.stepLengthProfile.plt");
  SteplengthProfile.output(o);
  o.close();
 
}

void TstVASteppingAction::UserSteppingAction(const G4Step* aStep)
{
  Steplength.accumulate(aStep->GetStepLength());

  G4double phi = aStep->GetDeltaPosition().phi();
  if (phi < 0.) phi = phi + twopi;
  SteplengthProfile.accumulate(phi,aStep->GetStepLength());

}
