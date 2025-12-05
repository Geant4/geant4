//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ActionInitializationMessenger.hh"
#include "ActionInitialization.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitializationMessenger::
ActionInitializationMessenger(ActionInitialization* act)
  :fAction(act)
{
  fActionDir = new G4UIdirectory("/action/");
  fActionDir->SetGuidance("Commands for action initialization status");

  fIAEAphspReaderDir = new G4UIdirectory("/action/IAEAphspReader/");
  fIAEAphspReaderDir->SetGuidance("Commands to set IAEAphsp reader object.");

  fIAEAphspWriterDir = new G4UIdirectory("/action/IAEAphspWriter/");
  fIAEAphspWriterDir->SetGuidance("Commands to set IAEAphsp writer object.");

  fIAEAphspReaderFileCmd =
    new G4UIcmdWithAString("/action/IAEAphspReader/fileName",this);
  fIAEAphspReaderFileCmd
    ->SetGuidance("Set IAEAphsp source file name, including path if needed.");
  fIAEAphspReaderFileCmd
    ->SetGuidance("(.IAEAphsp or .IAEAheader extension must not be written)");  
  fIAEAphspReaderFileCmd->SetParameterName("name",false);
  fIAEAphspReaderFileCmd->AvailableForStates(G4State_PreInit);

  fIAEAphspWriterFileCmd =
    new G4UIcmdWithAString("/action/IAEAphspWriter/namePrefix",this);
  fIAEAphspWriterFileCmd
    ->SetGuidance("Set the name prefix of IAEAphsp output files, ");
  fIAEAphspWriterFileCmd
    ->SetGuidance("including path if needed.");
  fIAEAphspWriterFileCmd
    ->SetGuidance("Name pattern: \"prefix_zphsp[_runID].IAEA[phsp|header]\".");
  fIAEAphspWriterFileCmd->SetGuidance("NOTE:");
  fIAEAphspWriterFileCmd
    ->SetGuidance("At least ONE zphsp value must be issued.");
  fIAEAphspWriterFileCmd->SetParameterName("prefix",false);
  fIAEAphspWriterFileCmd->AvailableForStates(G4State_PreInit);

  fIAEAphspWriterZphspCmd =
    new G4UIcmdWithADoubleAndUnit("/action/IAEAphspWriter/zphsp", this);
  fIAEAphspWriterZphspCmd
    ->SetGuidance("Add z-coordinate of output phsp plane.");
  fIAEAphspWriterZphspCmd
    ->SetGuidance("At least one value is NEEDED to get this writer working.");
  fIAEAphspWriterZphspCmd->SetParameterName("zphsp",false);
  fIAEAphspWriterZphspCmd->SetDefaultUnit("cm");
  fIAEAphspWriterZphspCmd->SetUnitCandidates("cm mm m");
  fIAEAphspWriterZphspCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ActionInitializationMessenger::~ActionInitializationMessenger()
{
  delete fActionDir;
  delete fIAEAphspReaderDir;
  delete fIAEAphspWriterDir;
  delete fIAEAphspReaderFileCmd;
  delete fIAEAphspWriterFileCmd;
  delete fIAEAphspWriterZphspCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ActionInitializationMessenger::SetNewValue(G4UIcommand* command,
						G4String newValue)
{
  if ( command == fIAEAphspReaderFileCmd )
    fAction->SetIAEAphspReader(newValue);

  else if ( command == fIAEAphspWriterFileCmd )
    fAction->SetIAEAphspWriterPrefix(newValue);

  else if ( command == fIAEAphspWriterZphspCmd )
    fAction->AddZphsp(fIAEAphspWriterZphspCmd->GetNewDoubleValue(newValue));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
