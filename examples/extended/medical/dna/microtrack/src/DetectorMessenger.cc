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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class

#include "DetectorMessenger.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* myDet)
  : fDetector(myDet)
{
  fDetectorDir = std::make_unique<G4UIdirectory>("/mygeom/");
  fDetectorDir->SetGuidance("Detector control.");

  fMaxRangeCmd =
    std::make_unique<G4UIcmdWithADoubleAndUnit>("/mygeom/maxRange", this);
  fMaxRangeCmd->SetGuidance("Maximum range of secondary electrons.");
  fMaxRangeCmd->SetGuidance("This is considered to set the z-size of ");
  fMaxRangeCmd->SetGuidance("the world volume.");
  fMaxRangeCmd->SetParameterName("range", false);
  fMaxRangeCmd->SetDefaultUnit("mm");
  fMaxRangeCmd->SetUnitCandidates("um mm cm");
  fMaxRangeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fHitSelRegZCmd =
    std::make_unique<G4UIcmdWithADoubleAndUnit>("/mygeom/hitSelRegZ", this);
  fHitSelRegZCmd->SetGuidance("Size of the hit selection region in Z");
  fHitSelRegZCmd->SetParameterName("side", false);
  fHitSelRegZCmd->SetDefaultUnit("mm");
  fHitSelRegZCmd->SetUnitCandidates("um mm cm");
  fHitSelRegZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fMatNameCmd = std::make_unique<G4UIcmdWithAString>("/mygeom/material", this);
  fMatNameCmd->SetGuidance("Name of the material.");
  fMatNameCmd->SetParameterName("name", false);
  fMatNameCmd->AvailableForStates(G4State_PreInit);

  fHitSelRegXYCmd =
    std::make_unique<G4UIcmdWithADoubleAndUnit>("/mygeom/hitSelRegXY", this);
  fHitSelRegXYCmd->SetGuidance("Size of the hit selection region in XY");
  fHitSelRegXYCmd->SetParameterName("side", false);
  fHitSelRegXYCmd->SetDefaultUnit("mm");
  fHitSelRegXYCmd->SetUnitCandidates("um mm cm");
  fHitSelRegXYCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fSiteRadiusCmd =
    std::make_unique<G4UIcmdWithADoubleAndUnit>("/mygeom/siteRadius", this);
  fSiteRadiusCmd->SetGuidance("Radius of the Site");
  fSiteRadiusCmd->SetParameterName("radius", false);
  fSiteRadiusCmd->SetDefaultUnit("um");
  fSiteRadiusCmd->SetUnitCandidates("nm um mm cm");
  fSiteRadiusCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fMaxRangeCmd.get()) {
    fDetector->SetMaxRange(fMaxRangeCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fMatNameCmd.get()) {
    fDetector->SetMaterial(newValue);
  }
  else if (command == fHitSelRegZCmd.get()) {
    fDetector->SetHitSelRegZ(fHitSelRegZCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fHitSelRegXYCmd.get()) {
    fDetector->SetHitSelRegXY(fHitSelRegXYCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fSiteRadiusCmd.get()) {
    fDetector->SetSiteRadius(fSiteRadiusCmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......