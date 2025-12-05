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
/// \file ExGflashMessenger.cc
/// \brief Implementation of the ExGflashMessenger class

#include "ExGflashMessenger.hh"

#include "ExGflashDetectorConstruction.hh"

#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIdirectory.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflashMessenger::ExGflashMessenger(ExGflashDetectorConstruction* Det) : fDetector(Det)
{
  fExGflashDir = new G4UIdirectory("/exgflash/");
  fExGflashDir->SetGuidance(" Gflash example commands.");

  fVerbose = new G4UIcmdWithAnInteger("/exgflash/verbose", this);
  fVerbose->SetGuidance("set exglash verbosity");
  fVerbose->SetGuidance("0- silent, 1 - init, 2 - run, 3 - event");
  fVerbose->SetParameterName("ver", false);
  fVerbose->SetRange("ver >= 0");
  fVerbose->AvailableForStates(G4State_PreInit, G4State_Idle);
  fVerbose->SetToBeBroadcasted(false);

  fDetDir = new G4UIdirectory("/exgflash/det/");
  fDetDir->SetGuidance("detector construction commands");

  fLBinCmd = new G4UIcmdWith3Vector("/exgflash/det/setLbin", this);
  fLBinCmd->SetGuidance("set longitudinal bining");
  fLBinCmd->SetGuidance("nLtot  - nb of bins; bin thickness (in Radl)");
  fLBinCmd->SetGuidance("dLradl - bin thickness (in fraction of Radl)");
  fLBinCmd->SetGuidance("LStart - dummy");
  fLBinCmd->SetParameterName("nLtot", "dLradl", "LStart", true);
  fLBinCmd->SetDefaultValue(G4ThreeVector{100.0, 0.25, 0.0});
  fLBinCmd->SetRange("nLtot>=1 && dLradl>0");
  fLBinCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fLBinCmd->SetToBeBroadcasted(false);

  fRBinCmd = new G4UIcmdWith3Vector("/exgflash/det/setRbin", this);
  fRBinCmd->SetGuidance("set radial bining");
  fRBinCmd->SetGuidance("nRtot  - nb of bins; bin thickness (in Rm)");
  fRBinCmd->SetGuidance("dRradl - bin thickness (in fraction of Rm)");
  fRBinCmd->SetGuidance("RStart - dummy");
  fRBinCmd->SetParameterName("nRtot", "dRradl", "RStart", true);
  fRBinCmd->SetDefaultValue(G4ThreeVector{80.0, 0.05, 0.0});
  fRBinCmd->SetRange("nRtot>=1 && dRradl>0");
  fRBinCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fRBinCmd->SetToBeBroadcasted(false);

  fMaterCmd = new G4UIcmdWithAString("/exgflash/det/setMat", this);
  fMaterCmd->SetGuidance("Select Material.");
  fMaterCmd->SetParameterName("material", false);
  fMaterCmd->AvailableForStates(G4State_PreInit);
  fMaterCmd->SetToBeBroadcasted(false);

  fNbCrysCmd = new G4UIcmdWithAnInteger("/exgflash/det/setNbCrys", this);
  fNbCrysCmd->SetGuidance("set numbber of crystals in row/col");
  fNbCrysCmd->SetGuidance("ncrys  - number of crystals ");
  fNbCrysCmd->SetParameterName("ncrys", false);
  fNbCrysCmd->SetRange("ncrys > 0");
  fNbCrysCmd->AvailableForStates(G4State_PreInit);
  fNbCrysCmd->SetToBeBroadcasted(false);

  fWCrysCmd = new G4UIcmdWithADoubleAndUnit("/exgflash/det/setCrysWidth", this);
  fWCrysCmd->SetGuidance("set crystal width (cm)");
  fWCrysCmd->SetGuidance("wcrys  - width of crystal ");
  fWCrysCmd->SetParameterName("wcrys", false);
  fWCrysCmd->SetRange("wcrys > 0");
  fWCrysCmd->SetUnitCategory("Length");
  fWCrysCmd->SetDefaultUnit("cm");
  fWCrysCmd->AvailableForStates(G4State_PreInit);
  fWCrysCmd->SetToBeBroadcasted(false);

  fLCrysCmd = new G4UIcmdWithADoubleAndUnit("/exgflash/det/setCrysLength", this);
  fLCrysCmd->SetGuidance("set crystal length (cm)");
  fLCrysCmd->SetGuidance("lcrys  - crystal length ");
  fLCrysCmd->SetParameterName("lcrys", false);
  fLCrysCmd->SetRange("lcrys > 0");
  fLCrysCmd->SetUnitCategory("Length");
  fLCrysCmd->SetDefaultUnit("cm");
  fLCrysCmd->AvailableForStates(G4State_PreInit);
  fLCrysCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExGflashMessenger::~ExGflashMessenger()
{
  delete fLCrysCmd;
  delete fWCrysCmd;
  delete fNbCrysCmd;
  delete fMaterCmd;
  delete fRBinCmd;
  delete fLBinCmd;
  delete fDetDir;
  delete fVerbose;
  delete fExGflashDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExGflashMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fMaterCmd) {
    fDetector->SetMaterial(newValue);
  }

  if (command == fNbCrysCmd) {
    fDetector->SetNbOfCrystals(fNbCrysCmd->GetNewIntValue(newValue));
  }

  if (command == fWCrysCmd) {
    fDetector->SetCrystalWidth(fWCrysCmd->GetNewDoubleValue(newValue));
  }

  if (command == fLCrysCmd) {
    fDetector->SetCrystalLength(fLCrysCmd->GetNewDoubleValue(newValue));
  }

  if (command == fVerbose) {
    fDetector->SetVerbose(fVerbose->GetNewIntValue(newValue));
  }

  if (command == fLBinCmd) {
    fDetector->SetLBining(fLBinCmd->GetNew3VectorValue(newValue));
  }

  if (command == fRBinCmd) {
    fDetector->SetRBining(fRBinCmd->GetNew3VectorValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
