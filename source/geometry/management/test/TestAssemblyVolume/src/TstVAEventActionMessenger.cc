// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: TstVAEventActionMessenger.cc,v 1.2 2001-02-01 21:27:14 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "TstVAEventActionMessenger.hh"

#include "TstVAEventAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TstVAEventActionMessenger::TstVAEventActionMessenger(TstVAEventAction* EvAct)
:eventAction(EvAct)
{ 
  DrawCmd = new G4UIcmdWithAString("/event/drawTracks",this);
  DrawCmd->SetGuidance("Draw the tracks in the event");
  DrawCmd->SetGuidance("  Choice : none, charged(default),neutral, all");
  DrawCmd->SetParameterName("choice",true);
  DrawCmd->SetDefaultValue("all");
  DrawCmd->SetCandidates("none charged neutral all");
  DrawCmd->AvailableForStates(Idle);
  
  PrintCmd = new G4UIcmdWithAnInteger("/event/printModulo",this);
  PrintCmd->SetGuidance("Print events modulo n");
  PrintCmd->SetParameterName("EventNb",false);
  PrintCmd->SetRange("EventNb>0");
  PrintCmd->AvailableForStates(Idle);     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TstVAEventActionMessenger::~TstVAEventActionMessenger()
{
  delete DrawCmd;
  delete PrintCmd;   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TstVAEventActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if(command == DrawCmd)
    {eventAction->SetDrawFlag(newValue);}
       
  if(command == PrintCmd)
    {eventAction->SetPrintModulo(PrintCmd->GetNewIntValue(newValue));}              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
