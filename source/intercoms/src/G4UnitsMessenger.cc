// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UnitsMessenger.cc,v 1.1 1999-01-07 16:09:28 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4UnitsMessenger.hh"

#include "G4UnitsTable.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4UnitsMessenger::G4UnitsMessenger()
{ 
  UnitsTableDir = new G4UIdirectory("/units/");
  UnitsTableDir->SetGuidance("Available units.");
      
  ListCmd = new G4UIcmdWithoutParameter("/units/list",this);
  ListCmd->SetGuidance("full list of available units.");
  ListCmd->AvailableForStates(Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4UnitsMessenger::~G4UnitsMessenger()
{
  delete ListCmd;
  delete UnitsTableDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4UnitsMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{  
  if (command == ListCmd)
    { G4UnitDefinition::PrintUnitsTable(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
