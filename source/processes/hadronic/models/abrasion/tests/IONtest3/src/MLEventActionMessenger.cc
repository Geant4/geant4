////////////////////////////////////////////////////////////////////////////////
//
#include "MLEventActionMessenger.hh"

#include "MLEventAction.hh"
////////////////////////////////////////////////////////////////////////////////
//
MLEventActionMessenger::MLEventActionMessenger (MLEventAction* EvAct)
  :eventAction(EvAct)
{ 
  DrawCmd = new G4UIcmdWithAString("/event/drawTracks",this);
  DrawCmd->SetGuidance("Draw the tracks in the event");
  DrawCmd->SetGuidance("Choice : none, charged(default), all");
  DrawCmd->SetParameterName("choice",true);
  DrawCmd->SetDefaultValue("charged");
  DrawCmd->SetCandidates("none charged all");
  DrawCmd->AvailableForStates(G4State_Idle);
  
  PrintCmd = new G4UIcmdWithAnInteger("/event/printModulo",this);
  PrintCmd->SetGuidance("Print events modulo n");
  PrintCmd->SetParameterName("EventNb",false);
  PrintCmd->SetRange("EventNb>0");
  PrintCmd->AvailableForStates(G4State_Idle);    

  CpuCmd = new G4UIcmdWithADouble("/run/cputime",this);
  CpuCmd->SetGuidance("Set the cputime limit (in seconds) for the run");
  CpuCmd->SetParameterName("cpu",false);
  CpuCmd->SetRange("cpu>0");
  CpuCmd->AvailableForStates(G4State_Idle);

  SeedCmd = new G4UIcmdWithAnInteger("/run/randseed",this);
  SeedCmd->SetGuidance(
    "Set the seed for random engine to one of the predefined values (0-214)");
  SeedCmd->SetParameterName("seed",false);
  SeedCmd->SetRange("seed>=0 && seed<215");
  SeedCmd->AvailableForStates(G4State_Idle);    
}
////////////////////////////////////////////////////////////////////////////////
//
//
MLEventActionMessenger::~MLEventActionMessenger ()
{
  delete DrawCmd;
  delete PrintCmd;  
  delete CpuCmd;
  delete SeedCmd;
}
////////////////////////////////////////////////////////////////////////////////
//
//
void MLEventActionMessenger::SetNewValue (G4UIcommand * command,
  G4String newValue)
{ 
  if (command == DrawCmd) {
    eventAction->SetDrawFlag(newValue);

  } else if (command == PrintCmd) {
    eventAction->SetPrintModulo(PrintCmd->GetNewIntValue(newValue));

  } else if (command == CpuCmd) {
    eventAction->SetCPUTime(CpuCmd->GetNewDoubleValue(newValue));

  } else if (command == SeedCmd) {
    eventAction->SetRandomSeed(SeedCmd->GetNewIntValue(newValue));

  }
}
////////////////////////////////////////////////////////////////////////////////
