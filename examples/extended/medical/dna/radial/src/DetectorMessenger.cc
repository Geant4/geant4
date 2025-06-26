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
/// \file DetectorMessenger.cc
/// \brief Implementation of the DetectorMessenger class

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* det, PhysicsList* pl)
:fDetector(det), fPhysList(pl)
{
  fDetDir = new G4UIdirectory("/radial/");
  fDetDir->SetGuidance("radial commands");

  fMatCmd = new G4UIcmdWithAString("/radial/setMat", this);
  fMatCmd->SetGuidance("Select material of the world.");
  fMatCmd->SetParameterName("Material", false);
  fMatCmd->AvailableForStates(G4State_PreInit);
  fMatCmd->SetToBeBroadcasted(false);

  fPhysCmd = new G4UIcmdWithAString("/radial/addPhysics", this);
  fPhysCmd->SetGuidance("Added Physics List");
  fPhysCmd->SetParameterName("Physics", false);
  fPhysCmd->AvailableForStates(G4State_PreInit);
  fPhysCmd->SetToBeBroadcasted(false);

  fTrackingCutCmd = new G4UIcmdWithABool("/radial/addIonsTrackingCut", this);
  fTrackingCutCmd->SetGuidance("Added Ions Tracking Cut");
  fTrackingCutCmd->SetDefaultValue(false);
  fTrackingCutCmd->AvailableForStates(G4State_PreInit);
  fTrackingCutCmd->SetToBeBroadcasted(false);

  fWorldRadiusCmd = new G4UIcmdWithADoubleAndUnit("/radial/setWorldRadius",this);
  fWorldRadiusCmd->SetGuidance("Set size of the World");
  fWorldRadiusCmd->SetParameterName("Size",false);
  fWorldRadiusCmd->SetRange("Size>0.");
  fWorldRadiusCmd->SetUnitCategory("Length");
  fWorldRadiusCmd->AvailableForStates(G4State_PreInit);

  fWorldLengthCmd = new G4UIcmdWithADoubleAndUnit("/radial/setWorldLength",this);
  fWorldLengthCmd->SetGuidance("Set size of the World");
  fWorldLengthCmd->SetParameterName("Size",false);
  fWorldLengthCmd->SetRange("Size>0.");
  fWorldLengthCmd->SetUnitCategory("Length");
  fWorldLengthCmd->AvailableForStates(G4State_PreInit);

  fThicknessCylindersCmd = new G4UIcmdWithADoubleAndUnit("/radial/setThicknessCylinders",this);
  fThicknessCylindersCmd->SetGuidance("Set thickness of cylinders");
  fThicknessCylindersCmd->SetParameterName("Size",false);
  fThicknessCylindersCmd->SetRange("Size>0.");
  fThicknessCylindersCmd->SetUnitCategory("Length");
  fThicknessCylindersCmd->AvailableForStates(G4State_PreInit);

  fMinRadiusCylindersCmd = new G4UIcmdWithADoubleAndUnit("/radial/setMinRadiusCylinders",this);
  fMinRadiusCylindersCmd->SetGuidance("Set minimum radius of the cylinders");
  fMinRadiusCylindersCmd->SetParameterName("Size",false);
  fMinRadiusCylindersCmd->SetRange("Size>0.");
  fMinRadiusCylindersCmd->SetUnitCategory("Length");
  fMinRadiusCylindersCmd->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fDetDir;

  delete fMatCmd;
  delete fPhysCmd;
  delete fTrackingCutCmd;
  delete fWorldRadiusCmd;
  delete fWorldLengthCmd;
  delete fThicknessCylindersCmd;
  delete fMinRadiusCylindersCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fMatCmd) fDetector->SetMaterial(newValue);

  if (command == fPhysCmd) fPhysList->AddPhysics(newValue);

  if (command == fTrackingCutCmd)
    fPhysList->SetTrackingCut(fTrackingCutCmd->GetNewBoolValue(newValue));

  if (command == fWorldRadiusCmd)
    fDetector->SetWorldRadius(fWorldRadiusCmd->GetNewDoubleValue(newValue));

  if (command == fWorldLengthCmd)
    fDetector->SetWorldLength(fWorldLengthCmd->GetNewDoubleValue(newValue));

  if (command == fThicknessCylindersCmd)
    fDetector->SetThicknessCylinders(fThicknessCylindersCmd->GetNewDoubleValue(newValue));

  if (command == fMinRadiusCylindersCmd)
    fDetector->SetMinRadiusCylinders(fMinRadiusCylindersCmd->GetNewDoubleValue(newValue));
}
