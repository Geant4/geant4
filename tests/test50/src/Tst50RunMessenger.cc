#include "Tst50RunMessenger.hh"
#include "Tst50RunAction.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"
Tst50RunMessenger::Tst50RunMessenger(Tst50RunAction* Tst50run)
:p_Run(Tst50run) 
{

  RunDir = new  G4UIdirectory("/run/test/");
  RunDir->SetGuidance("run GUI");
 


  testCmd= new G4UIcmdWithAString("/run/test/trans",this);
  testCmd->SetGuidance("chooce if you want transmission/backscattering test");
  testCmd->SetGuidance("Choice : on, off (default)");
  testCmd->SetParameterName("choice",true);
  testCmd->SetDefaultValue("off");
  testCmd->SetCandidates("on off");
  testCmd->AvailableForStates(G4State_Idle);


}
Tst50RunMessenger::~Tst50RunMessenger()
{
  delete  RunDir;
  delete  testCmd;
} 
void Tst50RunMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{

  if (command == testCmd){p_Run->Set_Trans(newValue);}

}
