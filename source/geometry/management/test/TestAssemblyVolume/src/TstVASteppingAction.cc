//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: TstVASteppingAction.cc,v 1.5 2003-06-16 16:52:18 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------

#include "TstVASteppingAction.hh"
#include "G4SteppingManager.hh"
#include <fstream>

TstVASteppingAction::TstVASteppingAction() : Steplength(100,0.,100.),
SteplengthProfile(100,0.,2*M_PI)
{ }

TstVASteppingAction::~TstVASteppingAction()
{
  std::ofstream o("test01.stepLength.plt");
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
