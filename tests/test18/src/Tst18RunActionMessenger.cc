// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		Tst18RunActionMessenger.cc
//
// Date:		16/08/99
// Author:		F Lei
// Organisation:	DERA UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		12115/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// 

#include "Tst18RunActionMessenger.hh"
#include "Tst18RunAction.hh"

#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst18RunActionMessenger::Tst18RunActionMessenger(Tst18RunAction* run)
:Tst18Run(run)
{ 
  FileCmd = new G4UIcmdWithAString("/run/filename",this);
  FileCmd->SetGuidance(" The log file name for the run");
  FileCmd->SetGuidance("   default = rdmex2.log ");
  FileCmd->SetParameterName(" Input ",true);
  FileCmd->SetDefaultValue("rdmex2.log");
  FileCmd->AvailableForStates(Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst18RunActionMessenger::~Tst18RunActionMessenger()
{
  delete FileCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst18RunActionMessenger::SetNewValue(G4UIcommand * command,G4String newValue)
{ 
    if(command == FileCmd)
    {Tst18Run->SetFilename(newValue);}
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....





