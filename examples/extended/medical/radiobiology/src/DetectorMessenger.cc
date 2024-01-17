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
/// \file radiobiology/src/DetectorMessenger.cc
/// \brief Implementation of the RadioBio::DetectorMessenger class

#include "DetectorMessenger.hh"

#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"

#include "DetectorConstruction.hh"

namespace RadioBio
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::DetectorMessenger(DetectorConstruction* Det) : fDetector(Det)
{
  fGeometryDir = new G4UIdirectory("/detectorGeom/");
  fGeometryDir->SetGuidance("commands to change geometry material and size");

  fMaterCmd = new G4UIcmdWithAString("/detectorGeom/setMat", this);
  fMaterCmd->SetGuidance("Select material of the box.");
  fMaterCmd->SetParameterName("choice", false);
  fMaterCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fMaterCmd->SetToBeBroadcasted(false);

  fSizeCmd = new G4UIcmdWithADoubleAndUnit("/detectorGeom/setSize", this);
  fSizeCmd->SetGuidance("Set size of the cubic box");
  fSizeCmd->SetParameterName("Size", false);
  fSizeCmd->SetRange("Size>0.");
  fSizeCmd->SetUnitCategory("Length");
  fSizeCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSizeCmd->SetToBeBroadcasted(false);

  fSizeVectorCmd = new G4UIcmdWith3VectorAndUnit("/detectorGeom/setBoxSizes", this);
  fSizeVectorCmd->SetGuidance("Insert sizes X Y and Z");
  fSizeVectorCmd->SetParameterName("SizeAlongX", "SizeAlongY", "SizeAlongZ", false);
  fSizeVectorCmd->SetRange("SizeAlongX>0. && SizeAlongY>0. && SizeAlongZ>0.");
  fSizeVectorCmd->SetUnitCategory("Length");
  fSizeVectorCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSizeVectorCmd->SetToBeBroadcasted(false);

  fSizeXCmd = new G4UIcmdWithADoubleAndUnit("/detectorGeom/setSizeX", this);
  fSizeXCmd->SetGuidance("Set X size of the box");
  fSizeXCmd->SetParameterName("Size", false);
  fSizeXCmd->SetRange("Size>0.");
  fSizeXCmd->SetUnitCategory("Length");
  fSizeXCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSizeXCmd->SetToBeBroadcasted(false);

  fSizeYCmd = new G4UIcmdWithADoubleAndUnit("/detectorGeom/setSizeY", this);
  fSizeYCmd->SetGuidance("Set Y size of the box");
  fSizeYCmd->SetParameterName("Size", false);
  fSizeYCmd->SetRange("Size>0.");
  fSizeYCmd->SetUnitCategory("Length");
  fSizeYCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSizeYCmd->SetToBeBroadcasted(false);

  fSizeZCmd = new G4UIcmdWithADoubleAndUnit("/detectorGeom/setSizeZ", this);
  fSizeZCmd->SetGuidance("Set Z size of the box");
  fSizeZCmd->SetParameterName("Size", false);
  fSizeZCmd->SetRange("Size>0.");
  fSizeZCmd->SetUnitCategory("Length");
  fSizeZCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fSizeZCmd->SetToBeBroadcasted(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorMessenger::~DetectorMessenger()
{
  delete fGeometryDir;
  delete fMaterCmd;
  delete fSizeCmd;
  delete fSizeVectorCmd;
  delete fSizeXCmd;
  delete fSizeYCmd;
  delete fSizeZCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fMaterCmd) {
    fDetector->SetMaterial(newValue);
  }

  if (command == fSizeCmd) {
    fDetector->SetSize(fSizeCmd->GetNewDoubleValue(newValue));
  }

  if (command == fSizeVectorCmd) {
    fDetector->SetSize(fSizeVectorCmd->GetNew3VectorValue(newValue));
  }

  if (command == fSizeXCmd) {
    fDetector->SetSizeX(fSizeXCmd->GetNewDoubleValue(newValue));
  }

  if (command == fSizeYCmd) {
    fDetector->SetSizeY(fSizeYCmd->GetNewDoubleValue(newValue));
  }

  if (command == fSizeZCmd) {
    fDetector->SetSizeZ(fSizeZCmd->GetNewDoubleValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace RadioBio