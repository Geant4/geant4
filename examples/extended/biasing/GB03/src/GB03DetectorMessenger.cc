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
/// \file GB03/src/GB03DetectorMessenger.cc
/// \brief Implementation of the GB03DetectorMessenger class

#include "GB03DetectorMessenger.hh"

#include "GB03DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4Material.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB03DetectorMessenger::GB03DetectorMessenger(GB03DetectorConstruction* det)
 : G4UImessenger(), fDetector(det)
{ 
  fDirectory = new G4UIdirectory("/GB03/");
  fDirectory->SetGuidance("UI commands of this example");
  
  G4String matList;
  const G4MaterialTable* matTbl = G4Material::GetMaterialTable();
  for(size_t i=0;i<G4Material::GetNumberOfMaterials();i++)
  {
    matList += (*matTbl)[i]->GetName();
    matList += " ";
  }

  fAbsMaterialCmd = new G4UIcmdWithAString("/GB03/setAbsMat",this);
  fAbsMaterialCmd->SetGuidance("Select Material of the Absorber.");
  fAbsMaterialCmd->SetParameterName("choice",false);
  fAbsMaterialCmd->AvailableForStates(G4State_Idle);
  fAbsMaterialCmd->SetCandidates(matList);
  
  fGapMaterialCmd = new G4UIcmdWithAString("/GB03/setGapMat",this);
  fGapMaterialCmd->SetGuidance("Select Material of the Gap.");
  fGapMaterialCmd->SetParameterName("choice",false);
  fGapMaterialCmd->AvailableForStates(G4State_Idle);
  fGapMaterialCmd->SetCandidates(matList);

  fNumLayerCmd = new G4UIcmdWithAnInteger("/GB03/numberOfLayers",this);
  fNumLayerCmd->SetGuidance("Set number of layers.");
  fNumLayerCmd->SetParameterName("nl",false);
  fNumLayerCmd->AvailableForStates(G4State_Idle);
  fNumLayerCmd->SetRange("nl>0");
    
  fVerboseCmd = new G4UIcmdWithAnInteger("/GB03/verbose",this);
  fVerboseCmd->SetGuidance("Set verbosity level");
  fVerboseCmd->SetParameterName("verbose",false);
  fVerboseCmd->AvailableForStates(G4State_Idle);
  fVerboseCmd->SetRange("verbose>=0");
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

GB03DetectorMessenger::~GB03DetectorMessenger()
{
  delete fAbsMaterialCmd;
  delete fGapMaterialCmd;
  delete fNumLayerCmd;
  delete fVerboseCmd;
  delete fDirectory;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void GB03DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fAbsMaterialCmd ) {
    fDetector->SetAbsorberMaterial(newValue);

  } else if( command == fGapMaterialCmd ){
    fDetector->SetGapMaterial(newValue);
  
  } else if( command == fNumLayerCmd ) {
    fDetector->SetNumberOfLayers(fNumLayerCmd->GetNewIntValue(newValue));

  } else if( command == fVerboseCmd ) {
    fDetector->SetVerboseLevel(fVerboseCmd->GetNewIntValue(newValue));
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String GB03DetectorMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String ans;
  if( command == fAbsMaterialCmd ){
    ans=fDetector->GetAbsorberMaterial(); 

  } else if( command == fGapMaterialCmd ){ 
    ans=fDetector->GetGapMaterial(); 

  } else if( command == fNumLayerCmd ) {
    ans=fNumLayerCmd->ConvertToString(GB03DetectorConstruction::GetNumberOfLayers()); 

  } else if( command == fVerboseCmd ) {
    ans=fVerboseCmd->ConvertToString(fDetector->GetVerboseLevel()); 
 
  }
  return ans;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
