// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em8EventActionMessenger.cc,v 1.1 2000-01-07 14:50:45 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em8EventActionMessenger.hh"

#include "Em8EventAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em8EventActionMessenger::Em8EventActionMessenger(Em8EventAction* EvAct)
:eventAction(EvAct)
{ 
  setVerboseCmd = new G4UIcmdWithAnInteger("/event/setverbose",this);
  setVerboseCmd->SetGuidance("Set verbose level ." );
  setVerboseCmd->SetParameterName("level",true);
  setVerboseCmd->SetDefaultValue(0);
  
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

Em8EventActionMessenger::~Em8EventActionMessenger()
{
  delete setVerboseCmd;
  delete DrawCmd;
  delete PrintCmd;    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em8EventActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if(command == setVerboseCmd)
    {eventAction->setEventVerbose(setVerboseCmd->GetNewIntValue(newValue));}
    
  if(command == DrawCmd)
    {eventAction->SetDrawFlag(newValue);}
    
  if(command == PrintCmd)
    {eventAction->SetPrintModulo(PrintCmd->GetNewIntValue(newValue));}               
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
