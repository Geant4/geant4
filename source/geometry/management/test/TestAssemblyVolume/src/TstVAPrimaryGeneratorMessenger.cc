
#include "TstVAPrimaryGeneratorMessenger.hh"

#include "TstVAPrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"
#include "G4ios.hh"

TstVAPrimaryGeneratorMessenger::TstVAPrimaryGeneratorMessenger(TstVAPrimaryGeneratorAction * myPGA)
:myPGAction(myPGA)
{
  G4bool omitable;

  G4UIdirectory* mydetDir = new G4UIdirectory("/myPGAction/");
  mydetDir->SetGuidance("Primary Generator Action control.");

  standardGun = new G4UIcmdWithoutParameter("/myPGAction/standardGun",this);
  standardGun->SetGuidance("Select standard gun.");

  randomGun = new G4UIcmdWithAString("/myPGAction/randomGun",this);
  randomGun->SetGuidance("Select random gun.");
  randomGun->SetParameterName("random-gun-type", omitable = true);
  randomGun->SetDefaultValue("randomDirection");
  randomGun->SetCandidates
    ("randomDirection randomPosition randomPositionAndDirection");
}

void TstVAPrimaryGeneratorMessenger::SetNewValue
(G4UIcommand * command,G4String newValues)
{
  if( command == standardGun)
  {
    myPGAction->SelectPrimaryGeneratorAction
      (TstVAPrimaryGeneratorAction::standardGun);
  }
  if( command == randomGun)
  {
    if (newValues == "randomDirection") {
      myPGAction->SelectPrimaryGeneratorAction
	(TstVAPrimaryGeneratorAction::randomDirectionGun);
    }
    if (newValues == "randomPosition") {
      myPGAction->SelectPrimaryGeneratorAction
	(TstVAPrimaryGeneratorAction::randomPositionGun);
    }
    if (newValues == "randomPositionAndDirection") {
      myPGAction->SelectPrimaryGeneratorAction
	(TstVAPrimaryGeneratorAction::randomPositionAndDirectionGun);
    }
  }
  return;
}
