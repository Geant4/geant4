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
/// \file src/DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorMessenger.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"

#include "DetectorConstruction.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det) : G4UImessenger(), fDetector(Det)
{
  fDetDir.reset(new G4UIdirectory("/det/"));
  fDetDir->SetGuidance("detector construction commands");

  fSizeYZCmd.reset(new G4UIcmdWithADoubleAndUnit("/det/setSizeYZ", this));
  fSizeYZCmd->SetGuidance("Set tranverse size of the calorimeter");
  fSizeYZCmd->SetParameterName("Size", false);
  fSizeYZCmd->SetRange("Size>0.");
  fSizeYZCmd->SetUnitCategory("Length");
  fSizeYZCmd->AvailableForStates(G4State_PreInit);
  fSizeYZCmd->SetToBeBroadcasted(false);

  fNbLayersCmd.reset(new G4UIcmdWithAnInteger("/det/setNbOfLayers", this));
  fNbLayersCmd->SetGuidance("Set number of layers.");
  fNbLayersCmd->SetParameterName("NbLayers", false);
  fNbLayersCmd->SetRange("NbLayers>0");
  fNbLayersCmd->AvailableForStates(G4State_PreInit);
  fNbLayersCmd->SetToBeBroadcasted(false);

  fNbAbsorCmd.reset(new G4UIcmdWithAnInteger("/det/setNbOfAbsor", this));
  fNbAbsorCmd->SetGuidance("Set number of Absorbers.");
  fNbAbsorCmd->SetParameterName("NbAbsor", false);
  fNbAbsorCmd->SetRange("NbAbsor>0");
  fNbAbsorCmd->AvailableForStates(G4State_PreInit);
  fNbAbsorCmd->SetToBeBroadcasted(false);

  fAbsorCmd.reset(new G4UIcommand("/det/setAbsor", this));
  fAbsorCmd->SetGuidance("Set the absor nb, the material, the thickness.");
  fAbsorCmd->SetGuidance("  absor number : from 1 to NbOfAbsor");
  fAbsorCmd->SetGuidance("  material name");
  fAbsorCmd->SetGuidance("  thickness (with unit) : t>0.");
  //
  G4UIparameter* AbsNbPrm = new G4UIparameter("AbsorNb", 'i', false);
  AbsNbPrm->SetGuidance("absor number : from 1 to NbOfAbsor");
  AbsNbPrm->SetParameterRange("AbsorNb>0");
  fAbsorCmd->SetParameter(AbsNbPrm);
  //
  G4UIparameter* MatPrm = new G4UIparameter("material", 's', false);
  MatPrm->SetGuidance("material name");
  fAbsorCmd->SetParameter(MatPrm);
  //
  G4UIparameter* ThickPrm = new G4UIparameter("thickness", 'd', false);
  ThickPrm->SetGuidance("thickness of absorber");
  ThickPrm->SetParameterRange("thickness>0.");
  fAbsorCmd->SetParameter(ThickPrm);
  //
  G4UIparameter* unitPrm = new G4UIparameter("unit", 's', false);
  unitPrm->SetGuidance("unit of thickness");
  G4String unitList = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("mm"));
  unitPrm->SetParameterCandidates(unitList);
  fAbsorCmd->SetParameter(unitPrm);
  //
  fAbsorCmd->AvailableForStates(G4State_PreInit);
  fAbsorCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fSizeYZCmd.get()) {
    fDetector->SetCalorSizeYZ(fSizeYZCmd->GetNewDoubleValue(newValue));
  }
  else if (command == fNbLayersCmd.get()) {
    fDetector->SetNbOfLayers(fNbLayersCmd->GetNewIntValue(newValue));
  }
  else if (command == fNbAbsorCmd.get()) {
    fDetector->SetNbOfAbsor(fNbAbsorCmd->GetNewIntValue(newValue));
  }
  else if (command == fAbsorCmd.get()) {
    G4int num;
    G4double tick;
    G4String unt, mat;
    std::istringstream is(newValue);
    is >> num >> mat >> tick >> unt;
    G4String material = mat;
    tick *= G4UIcommand::ValueOf(unt);
    fDetector->SetAbsorMaterial(num, material);
    fDetector->SetAbsorThickness(num, tick);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
