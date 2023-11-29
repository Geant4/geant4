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

// Author: Ivana Hrivnacova, 24/06/2013  (ivana@ipno.in2p3.fr)

#include "G4AnalysisMessenger.hh"
#include "G4VAnalysisManager.hh"
#include "G4NtupleMessenger.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4Threading.hh"

using namespace G4Analysis;

//_____________________________________________________________________________
G4AnalysisMessenger::G4AnalysisMessenger(G4VAnalysisManager* manager)
  : fManager(manager)
{
  fAnalysisDir = std::make_unique<G4UIdirectory>("/analysis/");
  fAnalysisDir->SetGuidance("analysis control");

  fOpenFileCmd = CreateCommand<G4UIcmdWithAString>(
    "openFile", "Open analysis file", "FileName", true);
  fOpenFileCmd->SetDefaultValue("");
  fOpenFileCmd->SetToBeBroadcasted(true);

  fWriteCmd = CreateCommandWithoutParameter(
    "write", "Write analysis data.");
  fWriteCmd->SetToBeBroadcasted(false);

  fResetCmd = CreateCommandWithoutParameter(
    "reset", "Reset analysis data.");
  fResetCmd->SetToBeBroadcasted(false);

  fCloseFileCmd = CreateCommand<G4UIcmdWithABool>(
    "closeFile", "Close analysis file and (optionally) reset data.", "IsReset", true);
  fCloseFileCmd->SetDefaultValue(true);
  fCloseFileCmd->SetToBeBroadcasted(false);

  fListCmd = CreateCommand<G4UIcmdWithABool>(
    "list", "List all/activate analysis objects.", "OnlyIfActive", true);
  fListCmd->SetDefaultValue(true);

  fSetActivationCmd = CreateCommand<G4UIcmdWithABool>(
    "setActivation",
    "Set activation. \n"
    "When this option is enabled, only the histograms marked as activated\n"
    "are returned, filled or saved on file.\n"
    "No warning is issued when Get or Fill is called on inactive histogram.",
    "Activation");

  fVerboseCmd = CreateCommand<G4UIcmdWithAnInteger>(
    "verbose", "Set verbose level", "VerboseLevel");
  fVerboseCmd->SetRange("VerboseLevel>=0 && VerboseLevel<=4");

  fCompressionCmd = CreateCommand<G4UIcmdWithAnInteger>(
    "compression", "Set compression level", "CompressionLevel");
  fCompressionCmd->SetRange("CompressionLevel>=0 && CompressionLevel<=4");

  fSetFileNameCmd = CreateCommand<G4UIcmdWithAString>(
    "setFileName", "Set name for the histograms & ntuple file", "Filename");

  fSetHistoDirNameCmd = CreateCommand<G4UIcmdWithAString>(
    "setHistoDirName", "Set name for the histograms directory", "HistoDirName");

  fSetNtupleDirNameCmd = CreateCommand<G4UIcmdWithAString>(
    "setNtupleDirName", "Set name for the ntuple directory", "NtupleDirName");

  fNtupleMessenger = std::make_unique<G4NtupleMessenger>(manager);
}

//_____________________________________________________________________________
G4AnalysisMessenger::~G4AnalysisMessenger() = default;

//
// private functions
//


//_____________________________________________________________________________
std::unique_ptr<G4UIcmdWithoutParameter>
G4AnalysisMessenger::CreateCommandWithoutParameter(
  G4String name, G4String guidance)
{
  G4String fullName = "/analysis/" + name;

  auto command = std::make_unique<G4UIcmdWithoutParameter>(fullName, this);
  command->SetGuidance(guidance.c_str());
  command->AvailableForStates(G4State_PreInit, G4State_Idle);

  return command;
}

//
// public functions
//

//_____________________________________________________________________________
void G4AnalysisMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if ( command == fOpenFileCmd.get() ) {
    // G4cout << "Calling OpenFile from command for " << fManager << G4endl;
    fManager->OpenFile(newValues);
    return;
  }

  if ( command == fWriteCmd.get() ) {
    fManager->WriteFromUI();
    return;
  }

  if ( command == fResetCmd.get() ) {
    fManager->ResetFromUI();
    return;
  }

  if ( command == fCloseFileCmd.get() ) {
    fManager->CloseFileFromUI(fCloseFileCmd->GetNewBoolValue(newValues));
    return;
  }

  if ( command == fListCmd.get() ) {
    fManager->List(fListCmd->GetNewBoolValue(newValues));
    return;
  }

  if ( command == fSetActivationCmd.get() ) {
    fManager->SetActivation(fSetActivationCmd->GetNewBoolValue(newValues));
    return;
  }

  if ( command == fVerboseCmd.get() ) {
    fManager->SetVerboseLevel(fVerboseCmd->GetNewIntValue(newValues));
    return;
  }

  if ( command == fCompressionCmd.get() ) {
    fManager->SetCompressionLevel(fCompressionCmd->GetNewIntValue(newValues));
    return;
  }

  if ( command == fSetFileNameCmd.get() ) {
    fManager->SetFileName(newValues);
    return;
  }

  if ( command == fSetHistoDirNameCmd.get() ) {
    fManager->SetHistoDirectoryName(newValues);
    return;
  }

  if ( command == fSetNtupleDirNameCmd.get() ) {
    fManager->SetNtupleDirectoryName(newValues);
    return;
  }
}
