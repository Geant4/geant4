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
// $Id$
//

#include "ExN07DetectorMessenger.hh"

#include "ExN07DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4Material.hh"

ExN07DetectorMessenger::ExN07DetectorMessenger(
                                           ExN07DetectorConstruction* ExN07Det)
:ExN07Detector(ExN07Det)
{ 
  N07Dir = new G4UIdirectory("/N07/");
  N07Dir->SetGuidance("UI commands of this example");
  
  G4String matList;
  const G4MaterialTable* matTbl = G4Material::GetMaterialTable();
  for(size_t i=0;i<G4Material::GetNumberOfMaterials();i++)
  {
    matList += (*matTbl)[i]->GetName();
    matList += " ";
  }

  AbsMaterCmd = new G4UIcmdWithAString("/N07/setAbsMat",this);
  AbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  AbsMaterCmd->SetParameterName("choice",false);
  AbsMaterCmd->AvailableForStates(G4State_Idle);
  AbsMaterCmd->SetCandidates(matList);
  
  GapMaterCmd = new G4UIcmdWithAString("/N07/setGapMat",this);
  GapMaterCmd->SetGuidance("Select Material of the Gap.");
  GapMaterCmd->SetParameterName("choice",false);
  GapMaterCmd->AvailableForStates(G4State_Idle);
  GapMaterCmd->SetCandidates(matList);

  numLayerCmd = new G4UIcmdWithAnInteger("/N07/numberOfLayers",this);
  numLayerCmd->SetGuidance("Set number of layers.");
  numLayerCmd->SetParameterName("nl",false);
  numLayerCmd->AvailableForStates(G4State_Idle);
  numLayerCmd->SetRange("nl>0");
    
  SerialCmd = new G4UIcmdWithABool("/N07/serialGeometry",this);
  SerialCmd->SetGuidance("Select calorimeters to be placed in serial or parallel.");
  SerialCmd->SetParameterName("serialize",false);
  SerialCmd->AvailableForStates(G4State_Idle);

  verboseCmd = new G4UIcmdWithAnInteger("/N07/verbose",this);
  verboseCmd->SetGuidance("Set verbosity level");
  verboseCmd->SetParameterName("verbose",false);
  verboseCmd->AvailableForStates(G4State_Idle);
  verboseCmd->SetRange("verbose>=0");
 
  AddMatCmd = new G4UIcmdWithABool("/N07/AddMaterial",this);
  AddMatCmd->SetGuidance("Add materials ");
  AddMatCmd->SetParameterName("dummy",true);
  AddMatCmd->AvailableForStates(G4State_Idle);


}

ExN07DetectorMessenger::~ExN07DetectorMessenger()
{
  delete AbsMaterCmd;
  delete GapMaterCmd;
  delete numLayerCmd;
  delete SerialCmd;
  delete N07Dir;  
}

void ExN07DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == AbsMaterCmd ) {
    ExN07Detector->SetAbsorberMaterial(newValue);

  } else if( command == GapMaterCmd ){
    ExN07Detector->SetGapMaterial(newValue);
  
  } else if( command == numLayerCmd ) {
    ExN07Detector->SetNumberOfLayers(numLayerCmd->GetNewIntValue(newValue));

  } else if( command == SerialCmd ) {
    ExN07Detector->SetSerialGeometry(SerialCmd->GetNewBoolValue(newValue));

  } else if( command == verboseCmd ) {
    ExN07Detector->SetVerboseLevel(verboseCmd->GetNewIntValue(newValue));

  } else if( command == AddMatCmd ) {
    ExN07Detector->AddMaterial();
    UpdateMaterialList(); 
  }
}

G4String ExN07DetectorMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String ans;
  if( command == AbsMaterCmd ){
    ans=ExN07Detector->GetAbsorberMaterial(); 

  } else if( command == GapMaterCmd ){ 
    ans=ExN07Detector->GetGapMaterial(); 

  } else if( command == numLayerCmd ) {
    ans=numLayerCmd->ConvertToString(ExN07Detector->GetNumberOfLayers()); 

  } else if( command == SerialCmd ){
    ans=SerialCmd->ConvertToString(ExN07Detector->IsSerial()); 

  } else if( command == SerialCmd ) {
    ans=SerialCmd->ConvertToString(ExN07Detector->IsSerial()); 
  
  } else if( command == verboseCmd ) {
    ans=verboseCmd->ConvertToString(ExN07Detector->GetVerboseLevel()); 
 
  }
  return ans;
}

void    ExN07DetectorMessenger::UpdateMaterialList() 
{
  G4String matList;
  const G4MaterialTable* matTbl = G4Material::GetMaterialTable();
  for(size_t i=0;i<G4Material::GetNumberOfMaterials();i++)
  {
    matList += (*matTbl)[i]->GetName();
    matList += " ";
  }

  if(AbsMaterCmd !=0) {
    AbsMaterCmd->SetCandidates(matList);
  }
  if (GapMaterCmd !=0) {
    GapMaterCmd->SetCandidates(matList);
  }
}


