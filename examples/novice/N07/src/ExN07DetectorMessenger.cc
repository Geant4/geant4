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
// $Id: ExN07DetectorMessenger.cc,v 1.2 2003/04/08 15:47:00 asaim Exp $
// GEANT4 tag $Name: geant4-05-01 $
//

#include "ExN07DetectorMessenger.hh"

#include "ExN07DetectorConstruction.hh"
#include "ExN07EventAction.hh"
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

  verboseCmd = new G4UIcmdWithAnInteger("/N07/verbose", this);
  verboseCmd->SetGuidance("Set verbosity of each event.");
  verboseCmd->SetParameterName("vl",true,true);
  verboseCmd->SetRange("vl>=0 && vl<10");
}

ExN07DetectorMessenger::~ExN07DetectorMessenger()
{
  delete AbsMaterCmd;
  delete GapMaterCmd;
  delete numLayerCmd;
  delete SerialCmd;
  delete verboseCmd;
  delete N07Dir;  
}

void ExN07DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == AbsMaterCmd )
  { ExN07Detector->SetAbsorberMaterial(newValue);}
  else if( command == GapMaterCmd )
  { ExN07Detector->SetGapMaterial(newValue);}
  else if( command == numLayerCmd )
  { ExN07Detector->SetNumberOfLayers(numLayerCmd->GetNewIntValue(newValue));}
  else if( command == SerialCmd )
  { ExN07Detector->SetSerialGeometry(SerialCmd->GetNewBoolValue(newValue));}
  else if( command == verboseCmd )
  { ExN07EventAction::SetVerboseLevel(verboseCmd->GetNewIntValue(newValue));}
}

G4String ExN07DetectorMessenger::GetCurrentValue(G4UIcommand * command)
{
  G4String ans;
  if( command == AbsMaterCmd )
  { ans=ExN07Detector->GetAbsorberMaterial(); }
  else if( command == GapMaterCmd )
  { ans=ExN07Detector->GetGapMaterial(); }
  else if( command == numLayerCmd )
  { ans=numLayerCmd->ConvertToString(ExN07Detector->GetNumberOfLayers()); }
  else if( command == SerialCmd )
  { ans=SerialCmd->ConvertToString(ExN07Detector->IsSerial()); }
  else if( command == verboseCmd )
  { ans=verboseCmd->ConvertToString(ExN07EventAction::GetVerboseLevel()); }
  return ans;
}



