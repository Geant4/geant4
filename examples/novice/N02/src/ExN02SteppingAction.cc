// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02SteppingAction.cc,v 1.2 1999-04-16 11:20:01 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "ExN02SteppingAction.hh"
#include "ExN02DetectorConstruction.hh"
#include "ExN02EventAction.hh"
#include "G4SteppingManager.hh"

ExN02SteppingAction::ExN02SteppingAction(ExN02DetectorConstruction* myDC, ExN02EventAction* myEA)
:myDetector(myDC), eventAction(myEA)
{ }

void ExN02SteppingAction::UserSteppingAction(const G4Step* aStep)
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


