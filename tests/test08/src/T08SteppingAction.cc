// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08SteppingAction.cc,v 1.2 1999-04-17 07:24:05 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "T08SteppingAction.hh"
#include "T08DetectorConstruction.hh"
#include "T08EventAction.hh"
#include "G4SteppingManager.hh"

T08SteppingAction::T08SteppingAction(T08DetectorConstruction* myDC, T08EventAction* myEA)
:myDetector(myDC), eventAction(myEA)
{ }

void T08SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // collect the energy deposited in the absorber
  
  const G4VPhysicalVolume* currentVolume = aStep->GetPreStepPoint()-> GetPhysicalVolume();

#if 0
  const G4VPhysicalVolume* absorber = myDetector->getAbsorber();

  if (currentVolume == absorber)
   {
    G4double EdepStep = aStep->GetTotalEnergyDeposit();
    eventAction->addEdep(EdepStep);
   } 
#endif
}


