// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MySteppingAction.cc,v 1.1 1999-01-07 16:05:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "MySteppingAction.hh"
#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "MySteppingActionMessenger.cc"

MySteppingAction::MySteppingAction()
:drawFlag(false)
{
  new MySteppingActionMessenger(this);
}

void MySteppingAction::UserSteppingAction()
{
  if(drawFlag)
  {
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

    if (pVVisManager) {
      const G4SteppingManager* pSM = GetSteppingManager();
      G4Polyline polyline;

      G4double charge = pSM->GetTrack()->GetDefinition()->GetPDGCharge();
      G4Colour colour;
      if      (charge < 0.) colour = G4Colour(1., 0., 0.);
      else if (charge > 0.) colour = G4Colour(0., 0., 1.);
      else                  colour = G4Colour(0., 1., 0.);
      G4VisAttributes attribs(colour);
      polyline.SetVisAttributes(attribs);
      polyline.append(pSM->GetStep()->GetPreStepPoint()->GetPosition());
      polyline.append(pSM->GetStep()->GetPostStepPoint()->GetPosition());
      pVVisManager -> Draw(polyline);
    }
  }
}


