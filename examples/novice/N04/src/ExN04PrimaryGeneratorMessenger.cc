
#include "ExN04PrimaryGeneratorMessenger.hh"
#include "ExN04PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"

ExN04PrimaryGeneratorMessenger::ExN04PrimaryGeneratorMessenger(ExN04PrimaryGeneratorAction * mpga)
:myAction(mpga)
{
  mydetDirectory = new G4UIdirectory("/mydet/");
  mydetDirectory->SetGuidance("ExN04 detector control commands.");

  genCmd = new G4UIcmdWithAString("/mydet/generator",this);
  genCmd->SetGuidance("Select primary generator.");
  genCmd->SetGuidance(" Available generators : PYTHIA, particleGun");
  genCmd->SetParameterName("generator",true);
  genCmd->SetDefaultValue("PYTHIA");
  genCmd->SetCandidates("PYTHIA particleGun");
}

ExN04PrimaryGeneratorMessenger::~ExN04PrimaryGeneratorMessenger()
{
  delete genCmd;
  delete mydetDirectory;
}

void ExN04PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{
  if( command==genCmd )
  { myAction->SetHEPEvtGenerator(newValue=="PYTHIA"); }
}

G4String ExN04PrimaryGeneratorMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String cv;
  
  if( command==genCmd )
  {
    if(myAction->GetHEPEvtGenerator())
    { cv = "PYTHIA"; }
    else
    { cv = "particleGun"; }
  }
  
  return cv;
}

