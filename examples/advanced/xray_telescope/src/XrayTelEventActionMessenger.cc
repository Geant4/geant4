//  XrayTelEventActionMessenger.cc
// 

#include "XrayTelEventActionMessenger.hh"

#include "XrayTelEventAction.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelEventActionMessenger::XrayTelEventActionMessenger(XrayTelEventAction* EvAct)
:eventAction(EvAct)
{ 
  DrawCmd = new G4UIcmdWithAString("/event/drawTracks",this);
  DrawCmd->SetGuidance("Draw the tracks in the event");
  DrawCmd->SetGuidance("  Choice : none, charged, all (default)");
  DrawCmd->SetParameterName("choice",true);
  DrawCmd->SetDefaultValue("charged");
  DrawCmd->SetCandidates("none charged all");
  DrawCmd->AvailableForStates(Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayTelEventActionMessenger::~XrayTelEventActionMessenger()
{
  delete DrawCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayTelEventActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
  if(command == DrawCmd)
    {eventAction->SetDrawFlag(newValue);}
}



		
