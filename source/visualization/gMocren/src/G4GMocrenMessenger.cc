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
// $Id: G4GMocrenMessenger.cc,v 1.1 2009-04-01 13:16:11 akimura Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Akinori Kimura    March 31, 2009
//
#include "G4GMocrenMessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"

G4GMocrenMessenger::G4GMocrenMessenger()
  : suffix (""), geometry(true), solids(true), invisibles(true),
    fgMocrenVolumeName("gMocrenVolume"),
    fgMocrenScoringMeshName("gMocrenScoringMesh"),
    fDrawVolumeGrid(false) {

  gMocrenDirectory = new G4UIdirectory("/vis/gMocren/");
  gMocrenDirectory->SetGuidance("gMocren commands.");

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

  fsetgMocrenVolumeNameCommand = new G4UIcmdWithAString("/vis/gMocren/setVolumeName", this);
  fsetgMocrenVolumeNameCommand->SetGuidance("detector name for a volume data in gMocren data.");
  fsetgMocrenVolumeNameCommand->SetParameterName("fgMocrenVolumeName",false);
  fsetgMocrenVolumeNameCommand->SetDefaultValue("gMocrenVolume");
  fsetgMocrenVolumeNameCommand->AvailableForStates(G4State_Idle);

  faddgMocrenHitNameCommand = new G4UIcmdWithAString("/vis/gMocren/addHitName", this);
  faddgMocrenHitNameCommand->SetGuidance("hit name for a dose distribution in gMocren data.");
  faddgMocrenHitNameCommand->SetParameterName("fgMocrenHitName",false);
  faddgMocrenHitNameCommand->AvailableForStates(G4State_Idle);

  fresetgMocrenHitNameCommand = new G4UIcmdWithoutParameter("/vis/gMocren/resetHitNames", this);
  fresetgMocrenHitNameCommand->SetGuidance("reset all hit names.");
  fresetgMocrenHitNameCommand->AvailableForStates(G4State_Idle);

  fsetgMocrenScoringMeshNameCommand = new G4UIcmdWithAString("/vis/gMocren/setScoringMeshName", this);
  fsetgMocrenScoringMeshNameCommand->SetGuidance("scoring mesh name for a dose distribution in gMocren data.");
  fsetgMocrenScoringMeshNameCommand->SetParameterName("fgMocrenScoringMeshName",false);
  fsetgMocrenScoringMeshNameCommand->SetDefaultValue("gMocrenScoringMesh");
  fsetgMocrenScoringMeshNameCommand->AvailableForStates(G4State_Idle);

  faddgMocrenHitScorerNameCommand = new G4UIcmdWithAString("/vis/gMocren/addHitScorerName", this);
  faddgMocrenHitScorerNameCommand->SetGuidance("hit scorer name for a dose distribution in gMocren data.");
  faddgMocrenHitScorerNameCommand->SetParameterName("fgMocrenHitScorerNames",false);
  faddgMocrenHitScorerNameCommand->AvailableForStates(G4State_Idle);

  fresetgMocrenHitScorerNameCommand = new G4UIcmdWithoutParameter("/vis/gMocren/resetHitScorerName", this);
  fresetgMocrenHitScorerNameCommand->SetGuidance("reset all hit scorer names.");
  fresetgMocrenHitScorerNameCommand->AvailableForStates(G4State_Idle);

  fsetgMocrenNoVoxelsCommand = new G4UIcommand("/vis/gMocren/setNumberOfVoxels", this);
  fsetgMocrenNoVoxelsCommand->SetGuidance("set number of voxels.");
  fsetgMocrenNoVoxelsCommand->AvailableForStates(G4State_Idle);
  G4UIparameter * param = new G4UIparameter("nX", 'i', false);
  param->SetDefaultValue("1");
  param->SetParameterRange("nX>0");
  fsetgMocrenNoVoxelsCommand->SetParameter(param);
  param = new G4UIparameter("nY", 'i', false);
  param->SetDefaultValue("1");
  param->SetParameterRange("nY>0");
  fsetgMocrenNoVoxelsCommand->SetParameter(param);
  param = new G4UIparameter("nZ", 'i', false);
  param->SetDefaultValue("1");
  param->SetParameterRange("nZ>0");
  fsetgMocrenNoVoxelsCommand->SetParameter(param);

  flistgMocrenCommand = new G4UIcmdWithoutParameter("/vis/gMocren/list", this);
  flistgMocrenCommand->SetGuidance("list gMocren command parameters.");
  flistgMocrenCommand->AvailableForStates(G4State_Idle);

  fDrawVolumeGridCommand = new G4UIcmdWithABool("/vis/gMocren/drawVolumeGrid", this);
  fDrawVolumeGridCommand->SetGuidance("Add grid of the volume.");
  fDrawVolumeGridCommand->SetParameterName("fDrawVolumeGrid",false);
  fDrawVolumeGridCommand->SetDefaultValue(false);
  fDrawVolumeGridCommand->AvailableForStates(G4State_Idle);

}

G4GMocrenMessenger::~G4GMocrenMessenger() {
  delete setEventNumberSuffixCommand;
  delete appendGeometryCommand;
  delete addPointAttributesCommand;
  delete useSolidsCommand;
  //    delete writeInvisiblesCommand;
  delete fsetgMocrenVolumeNameCommand;
  delete faddgMocrenHitNameCommand;
  delete fresetgMocrenHitNameCommand;
  //
  delete fsetgMocrenScoringMeshNameCommand;
  delete faddgMocrenHitScorerNameCommand;
  delete fresetgMocrenHitScorerNameCommand;
  //
  delete fsetgMocrenNoVoxelsCommand;
  //
  delete gMocrenDirectory;
  //
  delete fDrawVolumeGridCommand;
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
  } else if (command == fsetgMocrenVolumeNameCommand) {
    return fgMocrenVolumeName;
  } else if (command == faddgMocrenHitNameCommand) {
    G4String strval;
    std::vector<G4String>::iterator itr = fgMocrenHitNames.begin();
    for(; itr != fgMocrenHitNames.end(); itr++) {
      strval += *itr;
      strval += " ";
    }
    return strval;
  } else if (command == fsetgMocrenScoringMeshNameCommand) {
    return fgMocrenScoringMeshName;
  } else if (command == faddgMocrenHitScorerNameCommand) {
    G4String strval;
    std::vector<G4String>::iterator itr = fgMocrenHitScorerNames.begin();
    for(; itr != fgMocrenHitNames.end(); itr++) {
      strval += *itr;
      strval += " ";
    }
    return strval;
  } else if (command==fDrawVolumeGridCommand) {
    return fDrawVolumeGridCommand->ConvertToString(fDrawVolumeGrid);
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
  } else if (command == fsetgMocrenVolumeNameCommand) {
    fgMocrenVolumeName = newValue;
  } else if (command == faddgMocrenHitNameCommand) {
    fgMocrenHitNames.push_back(newValue);
  } else if (command == fresetgMocrenHitNameCommand) {
    fgMocrenHitNames.clear();
  } else if (command == fsetgMocrenScoringMeshNameCommand) {
    fgMocrenScoringMeshName = newValue;
  } else if (command == faddgMocrenHitScorerNameCommand) {
    fgMocrenHitScorerNames.push_back(newValue);
  } else if (command == fresetgMocrenHitScorerNameCommand) {
    fgMocrenHitScorerNames.clear();
  } else if (command == flistgMocrenCommand) {
    list();
  } else if (command == fsetgMocrenNoVoxelsCommand) {
    G4Tokenizer next(newValue);
    for(int i = 0; i < 3; i++) {
      fgMocrenNoVoxels[i] = StoI(next());
    }
  } else if (command==fDrawVolumeGridCommand) {
    fDrawVolumeGrid = fDrawVolumeGridCommand->GetNewBoolValue(newValue);
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
  return fgMocrenVolumeName;
}

std::vector<G4String> G4GMocrenMessenger::getHitNames() {
  return fgMocrenHitNames;
}

G4String G4GMocrenMessenger::getScoringMeshName() {
  return fgMocrenScoringMeshName;
}

std::vector<G4String> G4GMocrenMessenger::getHitScorerNames() {
  return fgMocrenHitScorerNames;
}

void G4GMocrenMessenger::list() {
  G4cout << "  Current valuess of gMocren command parameters:" << G4endl;
  //
  G4cout << "    volume name:        " << fgMocrenVolumeName << G4endl;
  //
  G4cout << "    hit names:          ";
  if(fgMocrenHitNames.size() > 0) {
    std::vector<G4String>::iterator itr = fgMocrenHitNames.begin();
    for(; itr != fgMocrenHitNames.end(); itr++)
      G4cout << *itr << "  " << G4endl;
  } else {
    G4cout << G4endl;
  }
  //
  G4cout << "    scoring mesh name:  " << fgMocrenScoringMeshName << G4endl;
  //
  G4cout << "    scorer names:       ";
  if(fgMocrenHitScorerNames.size() > 0) {
    std::vector<G4String>::iterator itr = fgMocrenHitScorerNames.begin();
    for(; itr != fgMocrenHitScorerNames.end(); itr++)
      G4cout << *itr << "  " << G4endl;
  } else {
    G4cout << G4endl;
  }
  G4cout << G4endl;
}

void G4GMocrenMessenger::getNoVoxels(G4int & nx, G4int & ny, G4int & nz) const {
  nx = fgMocrenNoVoxels[0];
  ny = fgMocrenNoVoxels[1];
  nz = fgMocrenNoVoxels[2];
}
