// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em0EventActionMessenger.cc,v 1.1 1999-01-08 16:32:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em0EventActionMessenger.hh"

#include "Em0EventAction.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em0EventActionMessenger::Em0EventActionMessenger(Em0EventAction* EvAct)
:eventAction(EvAct)
{ 
  DrawCmd = new G4UIcmdWithAString("/event/drawTracks",this);
  DrawCmd->SetGuidance("Draw the tracks in the event");
  DrawCmd->SetGuidance("  Choice : none,charged, all");
  DrawCmd->SetParameterName("choice",true);
  DrawCmd->SetDefaultValue("all");
  DrawCmd->SetCandidates("none charged all");
  DrawCmd->AvailableForStates(Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em0EventActionMessenger::~Em0EventActionMessenger()
{
  delete DrawCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em0EventActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if(command == DrawCmd)
    {eventAction->SetDrawFlag(newValue);}
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
