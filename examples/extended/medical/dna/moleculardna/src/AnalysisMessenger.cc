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
//
#include "AnalysisMessenger.hh"
#include "AnalysisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

AnalysisMessenger::AnalysisMessenger(AnalysisManager* analysisManager)
  : fpAnalysisManager(analysisManager)
  , fpAnalysisDirectory(new G4UIdirectory("/analysisDNA/"))
  , fpSaveStrands(new G4UIcmdWithABool("/analysisDNA/saveStrands", this))
  , fpStrandDirectory(new G4UIcmdWithAString("/analysisDNA/strandDir", this))
  , fpFragmentLength(new G4UIcmdWithAnInteger("/analysisDNA/fragmentGap", this))
  , fpSaveSingleChain(
      new G4UIcmdWithAnInteger("/analysisDNA/diagnosticChain", this))
  , fpDSBDistance(new G4UIcmdWithAnInteger("/analysisDNA/dsbDistance", this))
  , fpTestClassifier(
      new G4UIcmdWithoutParameter("/analysisDNA/testClassifier", this))
  , fpFileName(new G4UIcmdWithAString("/analysisDNA/fileName", this))
{
  // world geometry
  fpAnalysisDirectory->SetGuidance("App local commands for analysis");

  fpSaveStrands->SetGuidance("Bool to set whether strands ought be saved");
  fpSaveStrands->SetGuidance("use /analysisDNA/strandDir to set location");

  fpStrandDirectory->SetGuidance("Directory to save DNA damage fragments");
  fpStrandDirectory->SetParameterName("DNA framgents", false);

  fpFragmentLength->SetGuidance("Gap between DNA fragments in base pairs.");
  fpFragmentLength->SetGuidance(
    "Set to zero to score placement volumes independently");
  fpFragmentLength->SetParameterName("Base Pair gap", false);

  fpSaveSingleChain->SetGuidance(
    "Save the position of hits histos only on one chain");
  fpSaveSingleChain->SetParameterName("Chain Index", false);

  fpDSBDistance->SetGuidance("Max separation of DSBs. Must be less than 31.");
  fpDSBDistance->SetParameterName("Max. DSB distance.", false);

  fpTestClassifier->SetGuidance("Run unit test for break classification");

  fpFileName->SetGuidance("ROOT output file name");
  fpFileName->SetParameterName("ROOT output file name", false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void AnalysisMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fpStrandDirectory.get())
  {
    fpAnalysisManager->SetStrandDirectory(newValue);
  }
  else if(command == fpSaveStrands.get())
  {
    fpAnalysisManager->SetSaveStrands(
      G4UIcmdWithABool::GetNewBoolValue(newValue));
  }
  else if(command == fpFragmentLength.get())
  {
    fpAnalysisManager->SetFragmentGap(
      G4UIcmdWithAnInteger::GetNewIntValue(newValue));
  }
  else if(command == fpSaveSingleChain.get())
  {
    fpAnalysisManager->SetChainToSave(
      G4UIcmdWithAnInteger::GetNewIntValue(newValue));
  }
  else if(command == fpDSBDistance.get())
  {
    fpAnalysisManager->SetDSBDistance(
      G4UIcmdWithAnInteger::GetNewIntValue(newValue));
  }
  else if(command == fpTestClassifier.get())
  {
    fpAnalysisManager->TestClassification();
  }
  else if(command == fpFileName.get())
  {
    fpAnalysisManager->SetFileName(newValue);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
