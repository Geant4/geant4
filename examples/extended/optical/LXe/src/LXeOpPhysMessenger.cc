#include "LXeOpPhysMessenger.hh"
#include "LXeOpticalPhysics.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4Scintillation.hh"

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXeOpPhysMessenger::LXeOpPhysMessenger(LXeOpticalPhysics* LXeOpPhys)
:OpPhys(LXeOpPhys)
{
  yieldCmd = new G4UIcmdWithADouble("/LXe/detector/scintYieldFactor",this);
  yieldCmd->SetGuidance("Sets the yield factor for scintillation in all volumes");
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXeOpPhysMessenger::~LXeOpPhysMessenger()
{
  delete yieldCmd;
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXeOpPhysMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if (command == yieldCmd){
    OpPhys->SetScintYieldFactor(yieldCmd->GetNewDoubleValue(newValue));
  }
}


