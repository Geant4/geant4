#include "LXeSteppingMessenger.hh"
#include "LXeSteppingAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXeSteppingMessenger::LXeSteppingMessenger(LXeSteppingAction* step)
:stepping(step)
{
  oneStepPrimariesCmd = new G4UIcmdWithABool("/LXe/oneStepPrimaries",this);
  oneStepPrimariesCmd->SetGuidance("Only allows primaries to go one step in the scintillator volume before being killed.");
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXeSteppingMessenger::~LXeSteppingMessenger(){
  delete oneStepPrimariesCmd;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void 
LXeSteppingMessenger::SetNewValue(G4UIcommand* command,G4String newValue){ 
  if( command == oneStepPrimariesCmd ){ 
    stepping->SetOneStepPrimaries(oneStepPrimariesCmd
				  ->GetNewBoolValue(newValue));
  }
}


