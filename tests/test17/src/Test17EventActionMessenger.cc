// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Test17EventActionMessenger.hh"

#include "Test17EventAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17EventActionMessenger::Test17EventActionMessenger(Test17EventAction* EvAct)
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17EventActionMessenger::~Test17EventActionMessenger()
{
  delete setVerboseCmd;
  delete DrawCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17EventActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if(command == setVerboseCmd)
    {eventAction->setEventVerbose(setVerboseCmd->GetNewIntValue(newValue));}
    
  if(command == DrawCmd)
    {eventAction->SetDrawFlag(newValue);}    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
