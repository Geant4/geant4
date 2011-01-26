
#include "PhysicsListMessenger.hh"
#include "PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"


PhysicsListMessenger::PhysicsListMessenger(PhysicsList* physList) :
    physicsList(physList) {

  physicsDirectory = new G4UIdirectory("/physics/");
  physicsDirectory -> SetGuidance("Physics commands");

  physicsConstrCmd = 
               new G4UIcmdWithAString("/physics/physConstructor", this);  
  physicsConstrCmd -> SetGuidance("Registration of physics constructors");
  physicsConstrCmd -> SetParameterName("physConstructor",false);
  physicsConstrCmd -> AvailableForStates(G4State_PreInit);  

  prodThresholdCmd = 
               new G4UIcmdWithADoubleAndUnit("/physics/prodThreshold", this);
  prodThresholdCmd -> SetGuidance("Specification of production threshold");
  prodThresholdCmd -> SetParameterName("prodThreshold",true);
  prodThresholdCmd -> SetDefaultValue(0.001);
  prodThresholdCmd -> SetDefaultUnit("mm");
  prodThresholdCmd -> AvailableForStates(G4State_PreInit);  
}


PhysicsListMessenger::~PhysicsListMessenger() {

  delete prodThresholdCmd;
  delete physicsConstrCmd;
  delete physicsDirectory;
}


void PhysicsListMessenger::SetNewValue(G4UIcommand* cmd, G4String val) {

  if(cmd == physicsConstrCmd) 
     physicsList -> RegisterPhysConstructor(val);
  if(cmd == prodThresholdCmd) 
     physicsList -> SetProdThreshold(prodThresholdCmd -> GetNewDoubleValue(val));
}
