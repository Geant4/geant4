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
/// \file radiobiology/src/DoseMessenger.cc
/// \brief Implementation of the RadioBio::DoseMessenger class

#include "DoseMessenger.hh"

#include "Dose.hh"
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithADouble.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4UIdirectory.hh>

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DoseMessenger::DoseMessenger(Dose* dose) : G4UImessenger(), fDose(dose)
{
  // Diretory for dose commands
  fDoseDir = new G4UIdirectory("/dose/");
  fDoseDir->SetGuidance("commands to setup dose calculation");

  // Activate dose calculation
  fCalculationCmd = new G4UIcmdWithABool("/dose/calculate", this);
  fCalculationCmd->SetGuidance("Whether to enable dose calculation");
  fCalculationCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fCalculationCmd->SetToBeBroadcasted(false);

  // Set dose verbosity
  fVerbosityCmd = new G4UIcmdWithAnInteger("/dose/verbose", this);
  fVerbosityCmd->SetGuidance("Set verbosity level of dose");
  fVerbosityCmd->SetGuidance("0 = quiet");
  fVerbosityCmd->SetGuidance("1 = important messages (~10 per run)");
  fVerbosityCmd->SetGuidance("2 = debug");
  fVerbosityCmd->SetToBeBroadcasted(false);

  // Reset dose data
  fResetCmd = new G4UIcmdWithoutParameter("/dose/reset", this);
  fResetCmd->SetGuidance("Reset accumulated data");
  fResetCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fResetCmd->SetToBeBroadcasted(false);

  // Set dose filename
  fDosePathCmd = new G4UIcmdWithAString("/dose/fileName", this);
  fDosePathCmd->SetGuidance("Set the filename for the dose file");
  fDosePathCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fDosePathCmd->SetToBeBroadcasted(false);

  // Print Parameters
  fPrintCmd = new G4UIcmdWithoutParameter("/dose/print", this);
  fPrintCmd->SetGuidance("Print dose parameters");
  fPrintCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fPrintCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DoseMessenger::~DoseMessenger()
{
  delete fDoseDir;
  delete fCalculationCmd;
  delete fVerbosityCmd;
  delete fResetCmd;
  delete fDosePathCmd;
  delete fPrintCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DoseMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fCalculationCmd) {
    fDose->SetCalculationEnabled(fCalculationCmd->GetNewBoolValue(newValue));
  }

  if (command == fVerbosityCmd) {
    fDose->SetVerboseLevel(fVerbosityCmd->GetNewIntValue(newValue));
  }

  if (command == fResetCmd) {
    fDose->Reset();
  }

  if (command == fDosePathCmd) {
    fDose->SetPath(newValue);
  }

  if (command == fPrintCmd) {
    fDose->PrintParameters();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio