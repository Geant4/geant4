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
// $Id: G4GMocrenMessenger.cc 87126 2014-11-25 08:58:14Z gcosmo $
//
//
// Created:  Mar. 31, 2009  Akinori Kimura  
//
#include "G4GMocrenMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"

G4GMocrenMessenger::G4GMocrenMessenger()
  : suffix (""), geometry(true), pointAttributes(false), solids(true), invisibles(true),
    kgMocrenVolumeName("gMocrenVolume"),
    kgMocrenScoringMeshName("gMocrenScoringMesh"),
    kDrawVolumeGrid(false) {

  kgMocrenDirectory = new G4UIdirectory("/vis/gMocren/");
  kgMocrenDirectory->SetGuidance("gMocren commands.");

  setEventNumberSuffixCommand = new G4UIcmdWithAString("/vis/gMocren/setEventNumberSuffix", this);
  setEventNumberSuffixCommand->SetGuidance("Write separate event files, appended with given suffix.");
  setEventNumberSuffixCommand->SetGuidance("Define the suffix with a pattern such as '-0000'.");
  setEventNumberSuffixCommand->SetParameterName("suffix",false);
  setEventNumberSuffixCommand->SetDefaultValue("");
  setEventNumberSuffixCommand->AvailableForStates(G4State_Idle);
    
  appendGeometryCommand = new G4UIcmdWithABool("/vis/gMocren/appendGeometry", this);
  appendGeometryCommand->SetGuidance("Appends copy of geometry to every event.");
  appendGeometryCommand->SetParameterName("flag",false);
  appendGeometryCommand->SetDefaultValue(true);
  appendGeometryCommand->AvailableForStates(G4State_Idle);

  addPointAttributesCommand = new G4UIcmdWithABool("/vis/gMocren/addPointAttributes", this);
  addPointAttributesCommand->SetGuidance("Adds point attributes to the points of trajectories.");
  addPointAttributesCommand->SetParameterName("flag",false);
  addPointAttributesCommand->SetDefaultValue(false);
  addPointAttributesCommand->AvailableForStates(G4State_Idle);

  useSolidsCommand = new G4UIcmdWithABool("/vis/gMocren/useSolids", this);
  useSolidsCommand->SetGuidance("Use GMocren Solids, rather than Geant4 Primitives.");
  useSolidsCommand->SetParameterName("flag",false);
  useSolidsCommand->SetDefaultValue(true);
  useSolidsCommand->AvailableForStates(G4State_Idle);

  /* Not Enabled Yet
     writeInvisiblesCommand = new G4UIcmdWithABool("/vis/gMocren/writeInvisibles", this);
     writeInvisiblesCommand->SetGuidance("Write invisible objects.");
     writeInvisiblesCommand->SetParameterName("flag",false);
     writeInvisiblesCommand->SetDefaultValue(true);
     writeInvisiblesCommand->AvailableForStates(G4State_Idle);
  */

  kSetgMocrenVolumeNameCommand = new G4UIcmdWithAString("/vis/gMocren/setVolumeName", this);
  kSetgMocrenVolumeNameCommand->SetGuidance("detector name for a volume data in gMocren data.");
  kSetgMocrenVolumeNameCommand->SetParameterName("kgMocrenVolumeName",false);
  kSetgMocrenVolumeNameCommand->SetDefaultValue("gMocrenVolume");
  kSetgMocrenVolumeNameCommand->AvailableForStates(G4State_Idle);

  kAddgMocrenHitNameCommand = new G4UIcmdWithAString("/vis/gMocren/addHitName", this);
  kAddgMocrenHitNameCommand->SetGuidance("hit name for a dose distribution in gMocren data.");
  kAddgMocrenHitNameCommand->SetParameterName("kgMocrenHitName",false);
  kAddgMocrenHitNameCommand->AvailableForStates(G4State_Idle);

  kResetgMocrenHitNameCommand = new G4UIcmdWithoutParameter("/vis/gMocren/resetHitNames", this);
  kResetgMocrenHitNameCommand->SetGuidance("reset all hit names.");
  kResetgMocrenHitNameCommand->AvailableForStates(G4State_Idle);

  kSetgMocrenScoringMeshNameCommand = new G4UIcmdWithAString("/vis/gMocren/setScoringMeshName", this);
  kSetgMocrenScoringMeshNameCommand->SetGuidance("scoring mesh name for a dose distribution in gMocren data.");
  kSetgMocrenScoringMeshNameCommand->SetParameterName("kgMocrenScoringMeshName",false);
  kSetgMocrenScoringMeshNameCommand->SetDefaultValue("gMocrenScoringMesh");
  kSetgMocrenScoringMeshNameCommand->AvailableForStates(G4State_Idle);

  kAddgMocrenHitScorerNameCommand = new G4UIcmdWithAString("/vis/gMocren/addHitScorerName", this);
  kAddgMocrenHitScorerNameCommand->SetGuidance("hit scorer name for a dose distribution in gMocren data.");
  kAddgMocrenHitScorerNameCommand->SetParameterName("kgMocrenHitScorerNames",false);
  kAddgMocrenHitScorerNameCommand->AvailableForStates(G4State_Idle);

  kResetgMocrenHitScorerNameCommand = new G4UIcmdWithoutParameter("/vis/gMocren/resetHitScorerName", this);
  kResetgMocrenHitScorerNameCommand->SetGuidance("reset all hit scorer names.");
  kResetgMocrenHitScorerNameCommand->AvailableForStates(G4State_Idle);

  kSetgMocrenNoVoxelsCommand = new G4UIcommand("/vis/gMocren/setNumberOfVoxels", this);
  kSetgMocrenNoVoxelsCommand->SetGuidance("set number of voxels.");
  kSetgMocrenNoVoxelsCommand->AvailableForStates(G4State_Idle);
  G4UIparameter * param = new G4UIparameter("nX", 'i', false);
  param->SetDefaultValue("1");
  param->SetParameterRange("nX>0");
  kSetgMocrenNoVoxelsCommand->SetParameter(param);
  param = new G4UIparameter("nY", 'i', false);
  param->SetDefaultValue("1");
  param->SetParameterRange("nY>0");
  kSetgMocrenNoVoxelsCommand->SetParameter(param);
  param = new G4UIparameter("nZ", 'i', false);
  param->SetDefaultValue("1");
  param->SetParameterRange("nZ>0");
  kSetgMocrenNoVoxelsCommand->SetParameter(param);

  kListgMocrenCommand = new G4UIcmdWithoutParameter("/vis/gMocren/list", this);
  kListgMocrenCommand->SetGuidance("list gMocren command parameters.");
  kListgMocrenCommand->AvailableForStates(G4State_Idle);

  kDrawVolumeGridCommand = new G4UIcmdWithABool("/vis/gMocren/drawVolumeGrid", this);
  kDrawVolumeGridCommand->SetGuidance("Add grid of the volume.");
  kDrawVolumeGridCommand->SetParameterName("kDrawVolumeGrid",false);
  kDrawVolumeGridCommand->SetDefaultValue(false);
  kDrawVolumeGridCommand->AvailableForStates(G4State_Idle);

}

G4GMocrenMessenger::~G4GMocrenMessenger() {
  delete setEventNumberSuffixCommand;
  delete appendGeometryCommand;
  delete addPointAttributesCommand;
  delete useSolidsCommand;
  //    delete writeInvisiblesCommand;
  delete kSetgMocrenVolumeNameCommand;
  delete kAddgMocrenHitNameCommand;
  delete kResetgMocrenHitNameCommand;
  //
  delete kSetgMocrenScoringMeshNameCommand;
  delete kAddgMocrenHitScorerNameCommand;
  delete kResetgMocrenHitScorerNameCommand;
  //
  delete kSetgMocrenNoVoxelsCommand;
  //
  delete kgMocrenDirectory;
  //
  delete kDrawVolumeGridCommand;
}

G4String G4GMocrenMessenger::GetCurrentValue(G4UIcommand * command) {
  if (command==setEventNumberSuffixCommand) {
    return suffix;
  } else if (command==appendGeometryCommand) {
    return appendGeometryCommand->ConvertToString(geometry); 
  } else if (command==addPointAttributesCommand) {
    return addPointAttributesCommand->ConvertToString(pointAttributes); 
  } else if (command==useSolidsCommand) {
    return useSolidsCommand->ConvertToString(solids);
    //    } else if (command==writeInvisiblesCommand) {
    //        return writeInvisiblesCommand->ConvertToString(invisibles);
  } else if (command == kSetgMocrenVolumeNameCommand) {
    return kgMocrenVolumeName;
  } else if (command == kAddgMocrenHitNameCommand) {
    G4String strval;
    std::vector<G4String>::iterator itr = kgMocrenHitNames.begin();
    for(; itr != kgMocrenHitNames.end(); itr++) {
      strval += *itr;
      strval += " ";
    }
    return strval;
  } else if (command == kSetgMocrenScoringMeshNameCommand) {
    return kgMocrenScoringMeshName;
  } else if (command == kAddgMocrenHitScorerNameCommand) {
    G4String strval;
    std::vector<G4String>::iterator itr = kgMocrenHitScorerNames.begin();
    for(; itr != kgMocrenHitScorerNames.end(); itr++) {
      strval += *itr;
      strval += " ";
    }
    return strval;
  } else if (command==kDrawVolumeGridCommand) {
    return kDrawVolumeGridCommand->ConvertToString(kDrawVolumeGrid);
  } else {
    return "";
  }
}

void G4GMocrenMessenger::SetNewValue(G4UIcommand * command, G4String newValue) {
  if (command==setEventNumberSuffixCommand) {
    suffix = newValue;
  } else if (command==appendGeometryCommand) {
    geometry = appendGeometryCommand->GetNewBoolValue(newValue);
  } else if (command==addPointAttributesCommand) {
    pointAttributes = addPointAttributesCommand->GetNewBoolValue(newValue);
  } else if (command==useSolidsCommand) {
    solids = useSolidsCommand->GetNewBoolValue(newValue);
    //    } else if (command==writeInvisiblesCommand) {
    //        invisibles = writeInvisiblesCommand->GetNewBoolValue(newValue);
  } else if (command == kSetgMocrenVolumeNameCommand) {
    kgMocrenVolumeName = newValue;
  } else if (command == kAddgMocrenHitNameCommand) {
    kgMocrenHitNames.push_back(newValue);
  } else if (command == kResetgMocrenHitNameCommand) {
    kgMocrenHitNames.clear();
  } else if (command == kSetgMocrenScoringMeshNameCommand) {
    kgMocrenScoringMeshName = newValue;
  } else if (command == kAddgMocrenHitScorerNameCommand) {
    kgMocrenHitScorerNames.push_back(newValue);
  } else if (command == kResetgMocrenHitScorerNameCommand) {
    kgMocrenHitScorerNames.clear();
  } else if (command == kListgMocrenCommand) {
    list();
  } else if (command == kSetgMocrenNoVoxelsCommand) {
    G4Tokenizer next(newValue);
    for(int i = 0; i < 3; i++) {
      kgMocrenNoVoxels[i] = StoI(next());
    }
  } else if (command==kDrawVolumeGridCommand) {
    kDrawVolumeGrid = kDrawVolumeGridCommand->GetNewBoolValue(newValue);
  } 
}

G4String G4GMocrenMessenger::getEventNumberSuffix() {
  return suffix;
}

G4bool G4GMocrenMessenger::appendGeometry() {
  return geometry;
}

G4bool G4GMocrenMessenger::addPointAttributes() {
  return pointAttributes;
}

G4bool G4GMocrenMessenger::useSolids() {
  return solids;
}

G4bool G4GMocrenMessenger::writeInvisibles() {
  return invisibles;
}

G4String G4GMocrenMessenger::getVolumeName() {
  return kgMocrenVolumeName;
}

std::vector<G4String> G4GMocrenMessenger::getHitNames() {
  return kgMocrenHitNames;
}

G4String G4GMocrenMessenger::getScoringMeshName() {
  return kgMocrenScoringMeshName;
}

std::vector<G4String> G4GMocrenMessenger::getHitScorerNames() {
  return kgMocrenHitScorerNames;
}

void G4GMocrenMessenger::list() {
  G4cout << "  Current valuess of gMocren command parameters:" << G4endl;
  //
  G4cout << "    volume name:        " << kgMocrenVolumeName << G4endl;
  //
  G4cout << "    hit names:          ";
  if(kgMocrenHitNames.size() > 0) {
    std::vector<G4String>::iterator itr = kgMocrenHitNames.begin();
    for(; itr != kgMocrenHitNames.end(); itr++)
      G4cout << *itr << "  " << G4endl;
  } else {
    G4cout << G4endl;
  }
  //
  G4cout << "    scoring mesh name:  " << kgMocrenScoringMeshName << G4endl;
  //
  G4cout << "    scorer names:       ";
  if(kgMocrenHitScorerNames.size() > 0) {
    std::vector<G4String>::iterator itr = kgMocrenHitScorerNames.begin();
    for(; itr != kgMocrenHitScorerNames.end(); itr++)
      G4cout << *itr << "  " << G4endl;
  } else {
    G4cout << G4endl;
  }
  G4cout << G4endl;
}

void G4GMocrenMessenger::getNoVoxels(G4int & nx, G4int & ny, G4int & nz) const {
  nx = kgMocrenNoVoxels[0];
  ny = kgMocrenNoVoxels[1];
  nz = kgMocrenNoVoxels[2];
}
