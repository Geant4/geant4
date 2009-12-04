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

#include "ElectronEventAction.hh"


#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"

ElectronEventAction::ElectronEventAction()
{;}

ElectronEventAction::~ElectronEventAction()
{;}

void ElectronEventAction::BeginOfEventAction(const G4Event* event)
{
  if (fmod(event->GetEventID(), 100000.)==0)
	G4cout << "Event Number:" <<  event->GetEventID() << G4endl;
}

