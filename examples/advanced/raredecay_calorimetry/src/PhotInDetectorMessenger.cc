//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: PhotInDetectorMessenger.cc,v 1.2 2005-05-31 15:23:01 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#define debug

#include "PhotInDetectorMessenger.hh"

PhotInDetectorMessenger::PhotInDetectorMessenger(PhotInDetectorConstruction* PhotInDet):
PhotInDetector(PhotInDet)
{ 
#ifdef debug
  G4cout<<"PhotInDetectorMessenger::Constructor is called"<<G4endl;
#endif
  PhotInDir = new G4UIdirectory("/PhotIn/");
  PhotInDir->SetGuidance("UI commands of this example");
  
  AddMaterCmd = new G4UIcmdWithAString("/PhotIn/addMat",this);
  AddMaterCmd->SetGuidance("Add Material");
  AddMaterCmd->SetParameterName("choice",false);
  AddMaterCmd->AvailableForStates(G4State_PreInit);

  // Make a list of existing materials separated by a space
  G4String matList;
  const G4MaterialTable* matTbl = G4Material::GetMaterialTable();
  for(size_t i=0; i<G4Material::GetNumberOfMaterials(); i++)
  {
    matList += (*matTbl)[i]->GetName();
    matList += " ";
  }

  AbsMaterCmd = new G4UIcmdWithAString("/PhotIn/setAbsMat",this);
  AbsMaterCmd->SetGuidance("Select Material of the Absorber.");
  AbsMaterCmd->SetParameterName("choice",false);
  AbsMaterCmd->AvailableForStates(G4State_Idle);
  AbsMaterCmd->SetCandidates(matList);
  
  GapMaterCmd = new G4UIcmdWithAString("/PhotIn/setGapMat",this);
  GapMaterCmd->SetGuidance("Select Material of the Gap.");
  GapMaterCmd->SetParameterName("choice",false);
  GapMaterCmd->AvailableForStates(G4State_Idle);
  GapMaterCmd->SetCandidates(matList);

  numLayerCmd = new G4UIcmdWithAnInteger("/PhotIn/numberOfLayers",this);
  numLayerCmd->SetGuidance("Set number of layers.");
  numLayerCmd->SetParameterName("nl",false);
  numLayerCmd->AvailableForStates(G4State_Idle);
  numLayerCmd->SetRange("nl>0");

  numSlabsCmd = new G4UIcmdWithAnInteger("/PhotIn/numberOfSlabs",this);
  numSlabsCmd->SetGuidance("Set number of slabs in a layer.");
  numSlabsCmd->SetParameterName("ns",false);
  numSlabsCmd->AvailableForStates(G4State_Idle);
  numSlabsCmd->SetRange("ns>0");
    
  SerialCmd = new G4UIcmdWithABool("/PhotIn/serialGeometry",this);
  SerialCmd->SetGuidance("Select calorimeters to be placed in serial or parallel.");
  SerialCmd->SetParameterName("serialize",false);
  SerialCmd->AvailableForStates(G4State_Idle);

  verboseCmd = new G4UIcmdWithAnInteger("/PhotIn/verbose", this);
  verboseCmd->SetGuidance("Set verbosity of each event.");
  verboseCmd->SetParameterName("vl",true,true);
  verboseCmd->SetRange("vl>=0 && vl<10");
}

PhotInDetectorMessenger::~PhotInDetectorMessenger()
{
  delete AddMaterCmd;
  delete AbsMaterCmd;
  delete GapMaterCmd;
  delete numLayerCmd;
  delete numSlabsCmd;
  delete SerialCmd;
  delete verboseCmd;
  delete PhotInDir;  
}

void PhotInDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
#ifdef debug
  G4cout<<"PhotInDetectorMessenger::SetNewValue: "<<newValue<<G4endl;
#endif
  if( command == AddMaterCmd )
  {
    PhotInDetector->CreateMaterial(newValue);
    // update candidate list 
    G4String matList;
    const G4MaterialTable* matTbl = G4Material::GetMaterialTable();
    for(size_t i=0;i<G4Material::GetNumberOfMaterials();i++)
    {
      matList += (*matTbl)[i]->GetName();
      matList += " ";
    }
    GapMaterCmd->SetCandidates(matList);
    AbsMaterCmd->SetCandidates(matList);

  }
  else if( command == AbsMaterCmd ) PhotInDetector->SetAbsorberMaterial(newValue);
  else if( command == GapMaterCmd ) PhotInDetector->SetGapMaterial(newValue);
  else if( command == numLayerCmd )
    PhotInDetector->SetNumberOfLayers(numLayerCmd->GetNewIntValue(newValue));
  else if( command == numSlabsCmd )
    PhotInDetector->SetNumberOfSlabs(numSlabsCmd->GetNewIntValue(newValue));
  else if( command == SerialCmd )
    PhotInDetector->SetSerialGeometry(SerialCmd->GetNewBoolValue(newValue));
  else if( command == verboseCmd )
   PhotInEventAction::SetVerboseLevel(verboseCmd->GetNewIntValue(newValue));
  else G4cerr<<"***PhotInDetectorMessenger::SetNewValue: Command not found"<<G4endl;
}

G4String PhotInDetectorMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String ans;
  if( command == AbsMaterCmd ) ans=PhotInDetector->GetAbsorberMaterial();
  else if( command == GapMaterCmd ) ans=PhotInDetector->GetGapMaterial();
  else if( command == numLayerCmd )
    ans=numLayerCmd->ConvertToString(PhotInDetector->GetNumberOfLayers());
  else if( command == numSlabsCmd )
    ans=numSlabsCmd->ConvertToString(PhotInDetector->GetNumberOfSlabs());
  else if( command == SerialCmd )
    ans=SerialCmd->ConvertToString(PhotInDetector->IsSerial());
  else if( command == verboseCmd )
    ans=verboseCmd->ConvertToString(PhotInEventAction::GetVerboseLevel());
  else G4cerr<<"***PhotInDetectorMessenger::GetCurrentValue: Command not found"<<G4endl;
#ifdef debug
  G4cout<<"PhotInDetectorMessenger::GetCurrentValue: "<<ans<<G4endl;
#endif
  return ans;
}



