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
/// \file radiobiology/src/RBEMessenger.cc
/// \brief Implementation of the RadioBio::RBEMessenger class

#include "RBEMessenger.hh"

#include "RBE.hh"
#include <G4UIcmdWithABool.hh>
#include <G4UIcmdWithADouble.hh>
#include <G4UIcmdWithAString.hh>
#include <G4UIcmdWithAnInteger.hh>
#include <G4UIcmdWithoutParameter.hh>
#include <G4UIdirectory.hh>

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RBEMessenger::RBEMessenger(RBE* rbe) : G4UImessenger(), fRBE(rbe)
{
  // Directory for RBE commands
  fRBEDir = new G4UIdirectory("/rbe/");
  fRBEDir->SetGuidance("commands to setup RBE calculation");

  // Activate RBE calculation
  fCalculationCmd = new G4UIcmdWithABool("/rbe/calculate", this);
  fCalculationCmd->SetGuidance("Whether to enable RBE calculation");
  fCalculationCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fCalculationCmd->SetToBeBroadcasted(false);

  // Set RBE verbosity
  fVerbosityCmd = new G4UIcmdWithAnInteger("/rbe/verbose", this);
  fVerbosityCmd->SetGuidance("Set verbosity level of RBE");
  fVerbosityCmd->SetGuidance("0 = quiet");
  fVerbosityCmd->SetGuidance("1 = important messages (~10 per run)");
  fVerbosityCmd->SetGuidance("2 = debug");
  fVerbosityCmd->SetToBeBroadcasted(false);

  // Load a LEM table
  fLemTableCmd = new G4UIcmdWithAString("/rbe/loadLemTable", this);
  fLemTableCmd->SetGuidance("Load a LEM table used in calculating alpha&beta");
  fLemTableCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fLemTableCmd->SetToBeBroadcasted(false);

  // Load a cell line
  fCellLineCmd = new G4UIcmdWithAString("/rbe/cellLine", this);
  fCellLineCmd->SetGuidance("Set the cell line for alpha&beta calculation");
  fCellLineCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fCellLineCmd->SetToBeBroadcasted(false);

  // Reset RBE data
  fResetCmd = new G4UIcmdWithoutParameter("/rbe/reset", this);
  fResetCmd->SetGuidance("Reset accumulated data");
  fResetCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fResetCmd->SetToBeBroadcasted(false);

  // Print Parameters
  fPrintCmd = new G4UIcmdWithoutParameter("/rbe/print", this);
  fPrintCmd->SetGuidance("Print RBE parameters");
  fPrintCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fPrintCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RBEMessenger::~RBEMessenger()
{
  delete fRBEDir;
  delete fCalculationCmd;
  delete fVerbosityCmd;
  delete fLemTableCmd;
  delete fCellLineCmd;
  delete fResetCmd;
  delete fPrintCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RBEMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fCalculationCmd) {
    fRBE->SetCalculationEnabled(fCalculationCmd->GetNewBoolValue(newValue));
  }

  if (command == fVerbosityCmd) {
    fRBE->SetVerboseLevel(fVerbosityCmd->GetNewIntValue(newValue));
  }

  if (command == fLemTableCmd) {
    fRBE->LoadLEMTable(newValue);
  }

  if (command == fCellLineCmd) {
    fRBE->SetCellLine(newValue);
  }

  if (command == fResetCmd) {
    fRBE->Reset();
  }

  if (command == fPrintCmd) {
    fRBE->PrintParameters();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio