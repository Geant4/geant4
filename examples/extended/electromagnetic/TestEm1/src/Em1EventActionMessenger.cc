// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1EventActionMessenger.cc,v 1.1 1999-10-11 13:07:45 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em1EventActionMessenger.hh"

#include "Em1EventAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em1EventActionMessenger::Em1EventActionMessenger(Em1EventAction* EvAct)
:eventAction(EvAct)
{ 
  DrawCmd = new G4UIcmdWithAString("/event/drawTracks",this);
  DrawCmd->SetGuidance("Draw the tracks in the event");
  DrawCmd->SetGuidance("  Choice : none,charged, all");
  DrawCmd->SetParameterName("choice",true);
  DrawCmd->SetDefaultValue("all");
  DrawCmd->SetCandidates("none charged all");
  DrawCmd->AvailableForStates(Idle);
  
  PrintCmd = new G4UIcmdWithAnInteger("/event/printModulo",this);
  PrintCmd->SetGuidance("Print events modulo n");
  PrintCmd->SetParameterName("EventNb",false);
  PrintCmd->SetRange("EventNb>0");
  PrintCmd->AvailableForStates(Idle);      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em1EventActionMessenger::~Em1EventActionMessenger()
{
  delete DrawCmd;
  delete PrintCmd;   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1EventActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if(command == DrawCmd)
    {eventAction->SetDrawFlag(newValue);}
    
  if(command == PrintCmd)
    {eventAction->SetPrintModulo(PrintCmd->GetNewIntValue(newValue));}           
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
