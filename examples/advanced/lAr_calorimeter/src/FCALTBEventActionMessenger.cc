// FCALTBEventActionMessenger.cc v1.0
// /14/11/02 P.Mendez


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "FCALTBEventActionMessenger.hh"

#include "FCALTBEventAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FCALTBEventActionMessenger::FCALTBEventActionMessenger(FCALTBEventAction* EvAct)
:fcaltbeventAction(EvAct)
{ 
  DrawCmd = new G4UIcmdWithAString("/event/drawTracks",this);
  DrawCmd->SetGuidance("Draw the tracks in the event");
  DrawCmd->SetGuidance("  Choice : none,charged, all");
  DrawCmd->SetParameterName("choice",true);
  DrawCmd->SetDefaultValue("all");
  DrawCmd->SetCandidates("none charged all");
  DrawCmd->AvailableForStates(G4State_Idle);
  
  PrintCmd = new G4UIcmdWithAnInteger("/event/printModulo",this);
  PrintCmd->SetGuidance("Print events modulo n");
  PrintCmd->SetParameterName("EventNb",false);
  PrintCmd->SetRange("EventNb>0");
  PrintCmd->AvailableForStates(G4State_Idle);      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

FCALTBEventActionMessenger::~FCALTBEventActionMessenger()
{
  delete DrawCmd;
  delete PrintCmd;   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void FCALTBEventActionMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{ 
  if(command == DrawCmd)
    {fcaltbeventAction->SetDrawFlag(newValue);}
    
  if(command == PrintCmd)
    {fcaltbeventAction->SetPrintModulo(PrintCmd->GetNewIntValue(newValue));}           
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
