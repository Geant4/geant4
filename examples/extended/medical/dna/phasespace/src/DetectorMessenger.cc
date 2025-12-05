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

// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 51 (2024) 5873-5889
// Med. Phys. 45 (2018) e722-e739
// Phys. Med. 31 (2015) 861-874
// Med. Phys. 37 (2010) 4692-4708
// Int. J. Model. Simul. Sci. Comput. 1 (2010) 157–178
//
// The Geant4-DNA web site is available at http://geant4-dna.org
//

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det)
  : G4UImessenger(), fpDetector(Det)
{
  fpDetDir = new G4UIdirectory("/psp/");
  fpDetDir->SetGuidance("psp test commands");

  fpMaterCmd = new G4UIcmdWithAString("/psp/setMat", this);
  fpMaterCmd->SetGuidance("Select material of the world.");
  fpMaterCmd->SetParameterName("Material", false);
  fpMaterCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fpMaterCmd->SetToBeBroadcasted(false);

  fpDensityCmd = new G4UIcommand("/psp/setMatDens",this);
  fpDensityCmd->SetGuidance("Set density of the target material");
  G4UIparameter* symbPrm = new G4UIparameter("name",'s',false);
  symbPrm->SetGuidance("material name");
  fpDensityCmd->SetParameter(symbPrm);
  G4UIparameter* densityPrm = new G4UIparameter("density",'d',false);
  densityPrm->SetGuidance("density of material");
  densityPrm->SetParameterRange("density>0.");
  fpDensityCmd->SetParameter(densityPrm);
  G4UIparameter* unitPrm = new G4UIparameter("unit",'s',false);
  unitPrm->SetGuidance("unit of density");
  G4String unitList = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("g/cm3"));
  unitPrm->SetParameterCandidates(unitList);
  fpDensityCmd->SetParameter(unitPrm);
  fpDensityCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fpWorldSizeCmd = new G4UIcmdWithADoubleAndUnit("/psp/setWorldSize",this);
  fpWorldSizeCmd->SetGuidance("Set size of the World");
  fpWorldSizeCmd->SetParameterName("Size",false);
  fpWorldSizeCmd->SetRange("Size>0.");
  fpWorldSizeCmd->SetUnitCategory("Length");
  fpWorldSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fpScorerRadiusCmd = new G4UIcmdWithADoubleAndUnit("/psp/setScorerRadius",this);
  fpScorerRadiusCmd->SetGuidance("Set radius of the scoring sphere");
  fpScorerRadiusCmd->SetParameterName("Size",false);
  fpScorerRadiusCmd->SetRange("Size>0.");
  fpScorerRadiusCmd->SetUnitCategory("Length");
  fpScorerRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fpDetDir;

  delete fpMaterCmd;
  delete fpDensityCmd;
  delete fpWorldSizeCmd;
  delete fpScorerRadiusCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fpMaterCmd) fpDetector->SetMaterial(newValue);

  if (command == fpDensityCmd)
   {
     G4double dens;
     G4String name, unt;
     std::istringstream is(newValue);
     is >> name >> dens >> unt;
     dens *= G4UIcommand::ValueOf(unt);
     fpDetector->MaterialWithDensity(name,dens);
     fpDetector->SetMaterial(name);
   }

  if (command == fpWorldSizeCmd)
    fpDetector->SetWorldSize(fpWorldSizeCmd->GetNewDoubleValue(newValue));

  if (command == fpScorerRadiusCmd)
    fpDetector->SetScorerRadius(fpScorerRadiusCmd->GetNewDoubleValue(newValue));
}
