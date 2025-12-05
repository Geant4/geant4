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

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det, PhysicsList* PL)
  : G4UImessenger(), fpDetector(Det), fpPhysList(PL)
{
  fpDetDir = new G4UIdirectory("/dna/test/");
  fpDetDir->SetGuidance("dna test commands");

  fpMaterCmd = new G4UIcmdWithAString("/dna/test/setMat", this);
  fpMaterCmd->SetGuidance("Select material of the world.");
  fpMaterCmd->SetParameterName("Material", false);
  fpMaterCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fpMaterCmd->SetToBeBroadcasted(false);

  fpPhysCmd = new G4UIcmdWithAString("/dna/test/addPhysics", this);
  fpPhysCmd->SetGuidance("Added Physics List");
  fpPhysCmd->SetParameterName("Physics", false);
  fpPhysCmd->AvailableForStates(G4State_PreInit);
  fpPhysCmd->SetToBeBroadcasted(false);

  fpTrackingCutCmd = new G4UIcmdWithABool("/dna/test/addIonsTrackingCut", this);
  fpTrackingCutCmd->SetGuidance("Added Ions Tracking Cut");
  fpTrackingCutCmd->SetDefaultValue(false);
  fpTrackingCutCmd->AvailableForStates(G4State_PreInit);
  fpTrackingCutCmd->SetToBeBroadcasted(false);

  fDensityCmd = new G4UIcommand("/dna/test/setMatDens",this);
  fDensityCmd->SetGuidance("Set density of the target material");
  G4UIparameter* symbPrm = new G4UIparameter("name",'s',false);
  symbPrm->SetGuidance("material name");
  fDensityCmd->SetParameter(symbPrm);
  G4UIparameter* densityPrm = new G4UIparameter("density",'d',false);
  densityPrm->SetGuidance("density of material");
  densityPrm->SetParameterRange("density>0.");
  fDensityCmd->SetParameter(densityPrm);
  G4UIparameter* unitPrm = new G4UIparameter("unit",'s',false);
  unitPrm->SetGuidance("unit of density");
  G4String unitList = G4UIcommand::UnitsList(G4UIcommand::CategoryOf("g/cm3"));
  unitPrm->SetParameterCandidates(unitList);
  fDensityCmd->SetParameter(unitPrm);
  fDensityCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSizeCmd = new G4UIcmdWithADoubleAndUnit("/dna/test/setSize",this);
  fSizeCmd->SetGuidance("Set size of the World");
  fSizeCmd->SetParameterName("Size",false);
  fSizeCmd->SetRange("Size>0.");
  fSizeCmd->SetUnitCategory("Length");
  fSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fpDetDir;

  delete fpMaterCmd;
  delete fpPhysCmd;
  delete fpTrackingCutCmd;
  delete fDensityCmd;
  delete fSizeCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fpMaterCmd) fpDetector->SetMaterial(newValue);

  if (command == fpPhysCmd) fpPhysList->AddPhysics(newValue);

  if (command == fpTrackingCutCmd)
    fpPhysList->SetTrackingCut(fpTrackingCutCmd->GetNewBoolValue(newValue));

  if (command == fDensityCmd)
   {
     G4double dens;
     G4String name, unt;
     std::istringstream is(newValue);
     is >> name >> dens >> unt;
     dens *= G4UIcommand::ValueOf(unt);
     fpDetector->MaterialWithDensity(name,dens);
     fpDetector->SetMaterial(name);
   }

  if (command == fSizeCmd)
    fpDetector->SetSize(fSizeCmd->GetNewDoubleValue(newValue));
}
