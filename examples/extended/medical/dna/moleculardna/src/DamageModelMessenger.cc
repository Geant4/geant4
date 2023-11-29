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
#include "DamageModelMessenger.hh"

#include "DamageModel.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DamageModelMessenger::DamageModelMessenger(DamageModel* damageModel)
  : fpDamageModel(damageModel)
  , fpDamageDirectory(nullptr)
  , fpIntEnergyLower(nullptr)
  , fpIntEnergyUpper(nullptr)
  , fpOHBaseChance(nullptr)
  , fpOHStrandChance(nullptr)
  , fpOHInductionChance(nullptr)
  , fpHBaseChance(nullptr)
  , fpHStrandChance(nullptr)
  , fpHInductionChance(nullptr)
  , fpEaqBaseChance(nullptr)
  , fpEaqStrandChance(nullptr)
  , fpEaqInductionChance(nullptr)
{
  fpDamageDirectory = new G4UIdirectory("/dnadamage/");
  fpDamageDirectory->SetGuidance("Damage Model Parameters");

  fpIntEnergyLower =
    new G4UIcmdWithADoubleAndUnit("/dnadamage/directDamageLower", this);
  fpIntEnergyLower->SetGuidance("Minimum Energy required for an SSB");
  fpIntEnergyLower->SetGuidance(
    "Chance grows linearly until it reaches the upper value");
  fpIntEnergyLower->SetParameterName("Energy", false);

  fpIntEnergyUpper =
    new G4UIcmdWithADoubleAndUnit("/dnadamage/directDamageUpper", this);
  fpIntEnergyUpper->SetGuidance(
    "Energy required for an SSB to definitely occur");
  fpIntEnergyUpper->SetParameterName("Energy", false);

  // OH

  fpOHBaseChance =
    new G4UIcmdWithADouble("/dnadamage/indirectOHBaseChance", this);
  fpOHBaseChance->SetGuidance("Chance [0,1] of a OH damaging a base");
  fpOHBaseChance->SetParameterName("Energy", false);

  fpOHStrandChance =
    new G4UIcmdWithADouble("/dnadamage/indirectOHStrandChance", this);
  fpOHStrandChance->SetGuidance(
    "Chance [0,1] of a OH damaging sugar-phosphate");
  fpOHStrandChance->SetParameterName("Energy", false);

  fpOHInductionChance =
    new G4UIcmdWithADouble("/dnadamage/inductionOHChance", this);
  fpOHInductionChance->SetGuidance(
    "Chance [0,1] of a Base + OH -> strand break");
  fpOHInductionChance->SetParameterName("Energy", false);

  // H

  fpHBaseChance =
    new G4UIcmdWithADouble("/dnadamage/indirectHBaseChance", this);
  fpHBaseChance->SetGuidance("Chance [0,1] of a H damaging a base");
  fpHBaseChance->SetParameterName("Energy", false);

  fpHStrandChance =
    new G4UIcmdWithADouble("/dnadamage/indirectHStrandChance", this);
  fpHStrandChance->SetGuidance("Chance [0,1] of a H damaging sugar-phosphate");
  fpHStrandChance->SetParameterName("Energy", false);

  fpHInductionChance =
    new G4UIcmdWithADouble("/dnadamage/inductionHChance", this);
  fpHInductionChance->SetGuidance("Chance [0,1] of a Base + H -> strand break");
  fpHInductionChance->SetParameterName("Energy", false);

  // E_aq

  fpEaqBaseChance =
    new G4UIcmdWithADouble("/dnadamage/indirectEaqBaseChance", this);
  fpEaqBaseChance->SetGuidance("Chance [0,1] of a Eaq damaging a base");
  fpEaqBaseChance->SetParameterName("Energy", false);

  fpEaqStrandChance =
    new G4UIcmdWithADouble("/dnadamage/indirectEaqStrandChance", this);
  fpEaqStrandChance->SetGuidance(
    "Chance [0,1] of a Eaq damaging sugar-phosphate");
  fpEaqStrandChance->SetParameterName("Energy", false);

  fpEaqInductionChance =
    new G4UIcmdWithADouble("/dnadamage/inductionEaqChance", this);
  fpEaqInductionChance->SetGuidance(
    "Chance [0,1] of a Base + Eaq -> strand break");
  fpEaqInductionChance->SetParameterName("Energy", false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DamageModelMessenger::~DamageModelMessenger()
{
  // interaction model
  delete fpOHStrandChance;
  delete fpOHInductionChance;
  delete fpOHBaseChance;
  delete fpHStrandChance;
  delete fpHInductionChance;
  delete fpHBaseChance;
  delete fpEaqStrandChance;
  delete fpEaqInductionChance;
  delete fpEaqBaseChance;
  delete fpIntEnergyLower;
  delete fpIntEnergyUpper;
  delete fpIntRangeDirect;
  delete fpDamageDirectory;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DamageModelMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fpIntEnergyUpper)
  {
    fpDamageModel->SetDirectDamageUpper(
      G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue));
  }
  else if(command == fpIntEnergyLower)
  {
    fpDamageModel->SetDirectDamageLower(
      G4UIcmdWithADoubleAndUnit::GetNewDoubleValue(newValue));
  }
  else if(command == fpOHBaseChance)
  {
    fpDamageModel->SetIndirectOHBaseChance(
      G4UIcmdWithADouble::GetNewDoubleValue(newValue));
  }
  else if(command == fpOHStrandChance)
  {
    fpDamageModel->SetIndirectOHStrandChance(
      G4UIcmdWithADouble::GetNewDoubleValue(newValue));
  }
  else if(command == fpOHInductionChance)
  {
    fpDamageModel->SetInductionOHChance(
      G4UIcmdWithADouble::GetNewDoubleValue(newValue));
  }
  else if(command == fpHBaseChance)
  {
    fpDamageModel->SetIndirectHBaseChance(
      G4UIcmdWithADouble::GetNewDoubleValue(newValue));
  }
  else if(command == fpHStrandChance)
  {
    fpDamageModel->SetIndirectHStrandChance(
      G4UIcmdWithADouble::GetNewDoubleValue(newValue));
  }
  else if(command == fpHInductionChance)
  {
    fpDamageModel->SetInductionHChance(
      G4UIcmdWithADouble::GetNewDoubleValue(newValue));
  }
  else if(command == fpEaqBaseChance)
  {
    fpDamageModel->SetIndirectEaqBaseChance(
      G4UIcmdWithADouble::GetNewDoubleValue(newValue));
  }
  else if(command == fpEaqStrandChance)
  {
    fpDamageModel->SetIndirectEaqStrandChance(
      G4UIcmdWithADouble::GetNewDoubleValue(newValue));
  }
  else if(command == fpEaqInductionChance)
  {
    fpDamageModel->SetInductionEaqChance(
      G4UIcmdWithADouble::GetNewDoubleValue(newValue));
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
