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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunActionMessenger.hh"

#include "G4Tokenizer.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"

#include "RunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::RunActionMessenger(RunAction* runAction)
  : G4UImessenger(), fRunAction(runAction), fRunDir(nullptr), fScanParamsCmd(nullptr)
{
  fRunDir = new G4UIdirectory("/tomography/run/");
  fRunDir->SetGuidance("run control");

  fScanParamsCmd = new G4UIcommand("/tomography/run/scanParameters", this);
  fScanParamsCmd->SetGuidance("Set scan parameters: number of projections, slices and pixels");

  auto nbProjections = new G4UIparameter("NbProjections", 'i', false);
  nbProjections->SetGuidance("number of scan projections");
  nbProjections->SetParameterRange("NbProjections>=1");
  fScanParamsCmd->SetParameter(nbProjections);
  //
  auto nbSlices = new G4UIparameter("NbSlices", 'i', false);
  nbSlices->SetGuidance("number of scan slices");
  nbSlices->SetParameterRange("NbSlices>=1");
  fScanParamsCmd->SetParameter(nbSlices);
  //
  auto nbPixels = new G4UIparameter("NbPixels", 'i', false);
  nbPixels->SetGuidance("number of scan pixels");
  nbPixels->SetParameterRange("NbPixels>=1");
  fScanParamsCmd->SetParameter(nbPixels);
  fScanParamsCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  resumeCmd = new G4UIcmdWithABool("/tomography/run/resumeSimulation", this);
  resumeCmd->SetGuidance("Resume the simulation after an interruption");
  resumeCmd->SetParameterName("interruptionFlag", true);
  resumeCmd->SetDefaultValue(false);
  resumeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  resumeProjectionIndexCmd =
    new G4UIcmdWithAnInteger("/tomography/run/resumeProjectionIndex", this);
  resumeProjectionIndexCmd->SetGuidance("Set the index of projection to resume the simulation");
  resumeProjectionIndexCmd->SetParameterName("ProjectionIndex", true);
  resumeProjectionIndexCmd->SetDefaultValue(0);
  resumeProjectionIndexCmd->SetRange("ProjectionIndex >=0");
  resumeProjectionIndexCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunActionMessenger::~RunActionMessenger()
{
  delete fRunDir;
  delete fScanParamsCmd;
  delete resumeCmd;
  delete resumeProjectionIndexCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fScanParamsCmd) {
    G4Tokenizer next(newValue);
    G4int nbProjections = StoI(next());
    G4int nbSlices = StoI(next());
    G4int nbPixels = StoI(next());
    fRunAction->SetScanParameters(nbProjections, nbSlices, nbPixels);
  }
  else if (command == resumeCmd) {
    fRunAction->SetInterruptionFlag(resumeCmd->GetNewBoolValue(newValue));
  }
  else if (command == resumeProjectionIndexCmd) {
    fRunAction->SetResumeProjectionIndex(resumeProjectionIndexCmd->GetNewIntValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
