// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst01SteppingAction.cc,v 1.4 1999-12-15 14:54:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
////////////////////////////////////////////////////////////////////////////
//
//

#include "Tst01SteppingAction.hh"
#include "G4SteppingManager.hh"
#include "math.h"
#include "g4std/fstream"

////////////////////////////////////////////////////////////////////////////
//
//

Tst01SteppingAction::Tst01SteppingAction() : Steplength(100,0.,100.),
SteplengthProfile(100,0.,2*M_PI),fNumberOfTracks(10,0.,10.0)
{ }

//////////////////////////////////////////////////////////////////////////////
//
//

Tst01SteppingAction::~Tst01SteppingAction()
{
  G4std::ofstream o("NumberOfTracks.plt");
  fNumberOfTracks.output(o);
  o.close();


  //  G4std::ofstream o("test01.stepLength.plt");
  //  Steplength.output(o);
  //  o.close();
  //  o.open("test01.stepLengthProfile.plt");
  //  SteplengthProfile.output(o);
  //  o.close();
 
}



/////////////////////////////////////////////////////////////////////////
//
//

void Tst01SteppingAction::UserSteppingAction()
{
  G4Step* aStep = fpSteppingManager->GetStep() ;

  //  const G4Track* aTrack = GetSteppingManager()->GetTrack() ;

  // Get the number of steps ? from the track

  if( fpSteppingManager->GetTrack()->GetTrackStatus() != fAlive )
  {
     // Add this to an histogram of number of steps 

     fNumberOfTracks.accumulate((G4double)fpSteppingManager->GetTrack()->
                                                      GetCurrentStepNumber()) ;
  }

  // Steplength.accumulate(aStep->GetStepLength());

  // G4double phi = aStep->GetDeltaPosition().phi();
  // if (phi < 0.) phi = phi + twopi;
  // SteplengthProfile.accumulate(phi,aStep->GetStepLength());

}

//
//
/////////////////////////////////////////////////////////////////////////////





