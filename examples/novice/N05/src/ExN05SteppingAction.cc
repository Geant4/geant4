// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05SteppingAction.cc,v 1.3 1999-12-15 14:49:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "ExN05SteppingAction.hh"
#include "ExN05SteppingActionMessenger.hh"

#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SteppingManager.hh"

ExN05SteppingAction::ExN05SteppingAction()
:drawFlag(false)
{
  new ExN05SteppingActionMessenger(this);
}

void ExN05SteppingAction::UserSteppingAction(const G4Step*)
{
  if(drawFlag)
  {
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

    if (pVVisManager) {
      const G4SteppingManager* pSM = fpSteppingManager;
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


