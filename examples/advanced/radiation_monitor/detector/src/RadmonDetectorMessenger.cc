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
// File name:     RadmonDetectorMessenger.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMessenger.cc,v 1.3.2.2 2006/06/29 16:14:01 gunter Exp $
// Tag:           $Name: geant4-08-01 $
//

// Messenger commands path
#define COMMANDS_PATH "/radmon/detector/"

// Include files
#include "RadmonDetectorMessenger.hh"
#include "RadmonVDetectorLayout.hh"
#include "G4UnitsTable.hh"



                                                RadmonDetectorMessenger :: RadmonDetectorMessenger(RadmonVDetectorLayout * layout)
:
 RadmonMessenger(COMMANDS_PATH, "Interactive detector construction commands."),
 detectorLayout(layout),
 RADMON_INITIALIZE_COMMAND(EnableEnvironment),
 RADMON_INITIALIZE_COMMAND(DisableEnvironment),
 RADMON_INITIALIZE_COMMAND(SetEnvironmentType),
 RADMON_INITIALIZE_COMMAND(SetEnvironmentAttribute),
 RADMON_INITIALIZE_COMMAND(ClearEnvironmentAttribute),
 RADMON_INITIALIZE_COMMAND(CreateMultilayer),
 RADMON_INITIALIZE_COMMAND(RemoveMultilayer),
 RADMON_INITIALIZE_COMMAND(SetMultilayerWidth),
 RADMON_INITIALIZE_COMMAND(SetMultilayerHeight),
 RADMON_INITIALIZE_COMMAND(AppendLayerToMultilayer),
 RADMON_INITIALIZE_COMMAND(RemoveLayerFromMultilayer),
 RADMON_INITIALIZE_COMMAND(RemoveAllLayersFromMultilayer),
 RADMON_INITIALIZE_COMMAND(SetLayerThickness),
 RADMON_INITIALIZE_COMMAND(SetLayerType),
 RADMON_INITIALIZE_COMMAND(SetLayerAttribute),
 RADMON_INITIALIZE_COMMAND(ClearLayerAttribute),
 RADMON_INITIALIZE_COMMAND(CreatePlacement),
 RADMON_INITIALIZE_COMMAND(RemovePlacement),
 RADMON_INITIALIZE_COMMAND(SetPlacementPosition),
 RADMON_INITIALIZE_COMMAND(SetPlacementRotation),
 RADMON_INITIALIZE_COMMAND(SetRelativePlacementPosition),
 RADMON_INITIALIZE_COMMAND(SetRelativePlacementRotation),
 RADMON_INITIALIZE_COMMAND(DumpLayout),
 RADMON_INITIALIZE_COMMAND(Load),
 RADMON_INITIALIZE_COMMAND(Save)
{
 if (layout==0)
  G4Exception("RadmonDetectorMessenger::RadmonDetectorMessenger: layout==0.");

 RADMON_CREATE_COMMAND_0ARGS(EnableEnvironment,                 "Enables the environment");
 RADMON_CREATE_COMMAND_0ARGS(DisableEnvironment,                "Disables the environment");
 RADMON_CREATE_COMMAND_1ARG (SetEnvironmentType,                "Define the type of environment",                                                                       "type");
 RADMON_CREATE_COMMAND_2ARGS(SetEnvironmentAttribute,           "Define an attribute of the environment",                                                               "name", "value");
 RADMON_CREATE_COMMAND_1ARG (ClearEnvironmentAttribute,         "Removes an attribute of the environment",                                                              "name"); 
 RADMON_CREATE_COMMAND_1ARG (CreateMultilayer,                  "Creates a new multilayer",                                                                             "name");
 RADMON_CREATE_COMMAND_1ARG (RemoveMultilayer,                  "Removes a multilayer",                                                                                 "name");
 RADMON_CREATE_COMMAND_3ARGS(SetMultilayerWidth,                "Set the width of a multilayer",                                                                        "name", "width", "unit");
 RADMON_CREATE_COMMAND_3ARGS(SetMultilayerHeight,               "Set the height of a multilayer",                                                                       "name", "height", "unit");
 RADMON_CREATE_COMMAND_2ARGS(AppendLayerToMultilayer,           "Adds a layer to a multilayer",                                                                         "multilayerName", "layerName");
 RADMON_CREATE_COMMAND_2ARGS(RemoveLayerFromMultilayer,         "Removes a layer from a multilayer",                                                                    "multilayerName", "layerName");
 RADMON_CREATE_COMMAND_1ARG (RemoveAllLayersFromMultilayer,     "Removes all the layers of a multilayer",                                                               "multilayerName");
 RADMON_CREATE_COMMAND_4ARGS(SetLayerThickness,                 "Set the thickness of the layer of a multilayer",                                                       "multilayerName", "layerName", "width", "unit");
 RADMON_CREATE_COMMAND_3ARGS(SetLayerType,                      "Set the layer type of a multilayer",                                                                   "multilayerName", "layerName", "type");
 RADMON_CREATE_COMMAND_4ARGS(SetLayerAttribute,                 "Set an attribute of a layer of a multilayer",                                                          "multilayerName", "layerName", "attributeName", "value");
 RADMON_CREATE_COMMAND_3ARGS(ClearLayerAttribute,               "Removes an attribute of a layer of a multilayer",                                                      "multilayerName", "layerName", "attributeName");
 RADMON_CREATE_COMMAND_2ARGS(CreatePlacement,                   "Creates a new multilayer placement",                                                                   "placementName",  "multilayerName");
 RADMON_CREATE_COMMAND_1ARG (RemovePlacement,                   "Removes a multilayer placement",                                                                       "name");
 RADMON_CREATE_COMMAND_5ARGS(SetPlacementPosition,              "Changes the position of a placement",                                                                  "name", "x", "y", "z", "unit");
 RADMON_CREATE_COMMAND_5ARGS(SetPlacementRotation,              "Changes the orientation of a placement performing a rotation of delta along axis (theta, phi)",        "name", "theta", "phi", "delta", "unit");
 RADMON_CREATE_COMMAND_6ARGS(SetRelativePlacementPosition,      "Changes the position of a placement relative to another placement",                                    "name", "fromName", "x", "y", "z", "unit");
 RADMON_CREATE_COMMAND_6ARGS(SetRelativePlacementRotation,      "Changes the orientation of a placement relative to another placement",                                 "name", "fromName", "theta", "phi", "delta", "unit");
 RADMON_CREATE_COMMAND_0ARGS(DumpLayout,                        "Print out the current layout");
 RADMON_CREATE_COMMAND_1ARG (Load,                              "Loads a layout from file",                                                                             "fileName");
 RADMON_CREATE_COMMAND_1ARG (Save,                              "Saves a layout to file",                                                                               "fileName");
}



                                                RadmonDetectorMessenger :: ~RadmonDetectorMessenger()
{
 RADMON_DESTROY_COMMAND(Save);
 RADMON_DESTROY_COMMAND(Load);
 RADMON_DESTROY_COMMAND(DumpLayout);
 RADMON_DESTROY_COMMAND(SetRelativePlacementRotation);
 RADMON_DESTROY_COMMAND(SetRelativePlacementPosition);
 RADMON_DESTROY_COMMAND(SetPlacementRotation);
 RADMON_DESTROY_COMMAND(SetPlacementPosition);
 RADMON_DESTROY_COMMAND(RemovePlacement);
 RADMON_DESTROY_COMMAND(CreatePlacement);
 RADMON_DESTROY_COMMAND(ClearLayerAttribute);
 RADMON_DESTROY_COMMAND(SetLayerAttribute);
 RADMON_DESTROY_COMMAND(SetLayerType);
 RADMON_DESTROY_COMMAND(SetLayerThickness);
 RADMON_DESTROY_COMMAND(RemoveAllLayersFromMultilayer);
 RADMON_DESTROY_COMMAND(RemoveLayerFromMultilayer);
 RADMON_DESTROY_COMMAND(AppendLayerToMultilayer);
 RADMON_DESTROY_COMMAND(SetMultilayerHeight);
 RADMON_DESTROY_COMMAND(SetMultilayerWidth);
 RADMON_DESTROY_COMMAND(RemoveMultilayer);
 RADMON_DESTROY_COMMAND(CreateMultilayer);
 RADMON_DESTROY_COMMAND(ClearEnvironmentAttribute);
 RADMON_DESTROY_COMMAND(SetEnvironmentAttribute);
 RADMON_DESTROY_COMMAND(SetEnvironmentType);
 RADMON_DESTROY_COMMAND(DisableEnvironment);
 RADMON_DESTROY_COMMAND(EnableEnvironment);
}





G4String                                        RadmonDetectorMessenger :: GetCurrentValue(G4UIcommand * /* command */)
{
 G4cout << "RadmonDetectorMessenger::GetCurrentValue(): Not supported" << G4endl;
 
 return G4String();
}



void                                            RadmonDetectorMessenger :: SetNewValue(G4UIcommand * command, G4String newValue)
{
 RADMON_BEGIN_LIST_SET_COMMANDS
  RADMON_SET_COMMAND(EnableEnvironment)
  RADMON_SET_COMMAND(DisableEnvironment)
  RADMON_SET_COMMAND(SetEnvironmentType)
  RADMON_SET_COMMAND(SetEnvironmentAttribute)
  RADMON_SET_COMMAND(ClearEnvironmentAttribute)
  RADMON_SET_COMMAND(CreateMultilayer)
  RADMON_SET_COMMAND(RemoveMultilayer)
  RADMON_SET_COMMAND(SetMultilayerWidth)
  RADMON_SET_COMMAND(SetMultilayerHeight)
  RADMON_SET_COMMAND(AppendLayerToMultilayer)
  RADMON_SET_COMMAND(RemoveLayerFromMultilayer)
  RADMON_SET_COMMAND(RemoveAllLayersFromMultilayer)
  RADMON_SET_COMMAND(SetLayerThickness)
  RADMON_SET_COMMAND(SetLayerType)
  RADMON_SET_COMMAND(SetLayerAttribute)
  RADMON_SET_COMMAND(ClearLayerAttribute)
  RADMON_SET_COMMAND(CreatePlacement)
  RADMON_SET_COMMAND(RemovePlacement)
  RADMON_SET_COMMAND(SetPlacementPosition)
  RADMON_SET_COMMAND(SetPlacementRotation)
  RADMON_SET_COMMAND(SetRelativePlacementPosition)
  RADMON_SET_COMMAND(SetRelativePlacementRotation)
  RADMON_SET_COMMAND(DumpLayout)
  RADMON_SET_COMMAND(Load)
  RADMON_SET_COMMAND(Save)
 RADMON_END_LIST_SET_COMMANDS
}





// Events
void                                            RadmonDetectorMessenger :: OnEnableEnvironment(const G4String & /* value */)
{
 detectorLayout->EnableEnvironment();
}



void                                            RadmonDetectorMessenger :: OnDisableEnvironment(const G4String & /* value */)
{
 detectorLayout->DisableEnvironment();
}



void                                            RadmonDetectorMessenger :: OnSetEnvironmentType(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 

 detectorLayout->SetEnvironmentType(args[0]);
}



void                                            RadmonDetectorMessenger :: OnSetEnvironmentAttribute(const G4String & value)
{
 G4String args[2];

 if (!ProcessArguments(value, 2, args))
  return; 

 detectorLayout->SetEnvironmentAttribute(args[0], args[1]);
}



void                                            RadmonDetectorMessenger :: OnClearEnvironmentAttribute(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 

 detectorLayout->ClearEnvironmentAttribute(args[0]);
}



void                                            RadmonDetectorMessenger :: OnCreateMultilayer(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 

 detectorLayout->CreateMultilayer(args[0]);
}



void                                            RadmonDetectorMessenger :: OnRemoveMultilayer(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 

 detectorLayout->RemoveMultilayer(args[0]);
}



void                                            RadmonDetectorMessenger :: OnSetMultilayerWidth(const G4String & value)
{
 G4String args[3];

 if (!ProcessArguments(value, 3, args))
  return; 
  
 G4double dbl(GetUnit(args[2], "Length"));
 
 if (dbl<=0.)
  return;
 
 dbl*=G4UIcommand::ConvertToDouble(args[1]);

 if (dbl<0.)
 {
  G4cout << "RadmonDetectorMessenger::OnSetMultilayerWidth(): Cannot set negative width." << G4endl;
  return;
 }
 
 detectorLayout->SetMultilayerWidth(args[0], dbl);
}



void                                            RadmonDetectorMessenger :: OnSetMultilayerHeight(const G4String & value)
{
 G4String args[3];

 if (!ProcessArguments(value, 3, args))
  return; 
  
 G4double dbl(GetUnit(args[2], "Length"));
 
 if (dbl<=0.)
  return;
 
 dbl*=G4UIcommand::ConvertToDouble(args[1]);
 
 if (dbl<0)
 {
  G4cout << "RadmonDetectorMessenger::OnSetMultilayerHeight(): Cannot set negative height." << G4endl;
  return;
 }
 
 detectorLayout->SetMultilayerHeight(args[0], dbl);
}



void                                            RadmonDetectorMessenger :: OnAppendLayerToMultilayer(const G4String & value)
{
 G4String args[2];

 if (!ProcessArguments(value, 2, args))
  return; 

 detectorLayout->AppendLayerToMultilayer(args[0], args[1]);
}



void                                            RadmonDetectorMessenger :: OnRemoveLayerFromMultilayer(const G4String & value)
{
 G4String args[2];

 if (!ProcessArguments(value, 2, args))
  return; 

 detectorLayout->RemoveLayerFromMultilayer(args[0], args[1]);
}



void                                            RadmonDetectorMessenger :: OnRemoveAllLayersFromMultilayer(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 

 detectorLayout->RemoveAllLayersFromMultilayer(args[0]);
}



void                                            RadmonDetectorMessenger :: OnSetLayerThickness(const G4String & value)
{
 G4String args[4];

 if (!ProcessArguments(value, 4, args))
  return; 
  
 G4double dbl(GetUnit(args[3], "Length"));
 
 if (dbl<=0.)
  return;
 
 if (dbl<0.)
 {
  G4cout << "RadmonDetectorMessenger::OnSetLayerThickness(): Cannot set negative thickness." << G4endl;
  return;
 }
 
 dbl*=G4UIcommand::ConvertToDouble(args[2]);
 detectorLayout->SetLayerThickness(args[0], args[1], dbl);
}



void                                            RadmonDetectorMessenger :: OnSetLayerType(const G4String & value)
{
 G4String args[3];

 if (!ProcessArguments(value, 3, args))
  return; 

 detectorLayout->SetLayerType(args[0], args[1], args[2]);
}



void                                            RadmonDetectorMessenger :: OnSetLayerAttribute(const G4String & value)
{
 G4String args[4];

 if (!ProcessArguments(value, 4, args))
  return; 

 detectorLayout->SetLayerAttribute(args[0], args[1], args[2], args[3]);
}



void                                            RadmonDetectorMessenger :: OnClearLayerAttribute(const G4String & value)
{
 G4String args[3];

 if (!ProcessArguments(value, 3, args))
  return; 

 detectorLayout->ClearLayerAttribute(args[0], args[1], args[2]);
}



void                                            RadmonDetectorMessenger :: OnCreatePlacement(const G4String & value)
{
 G4String args[2];

 if (!ProcessArguments(value, 2, args))
  return; 

 detectorLayout->CreatePlacement(args[0], args[1]);
}



void                                            RadmonDetectorMessenger :: OnRemovePlacement(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 

 detectorLayout->RemovePlacement(args[0]);
}



void                                            RadmonDetectorMessenger :: OnSetPlacementPosition(const G4String & value)
{
 G4String args[5];

 if (!ProcessArguments(value, 5, args))
  return; 

 G4double x(GetUnit(args[4], "Length"));
 
 if (x<=0.)
  return;
 
 G4double y(x*G4UIcommand::ConvertToDouble(args[2]));
 G4double z(x*G4UIcommand::ConvertToDouble(args[3]));
 x*=G4UIcommand::ConvertToDouble(args[1]);

 detectorLayout->SetPlacementPosition(args[0], G4ThreeVector(x, y, z));
}



void                                            RadmonDetectorMessenger :: OnSetPlacementRotation(const G4String & value)
{
 G4String args[5];

 if (!ProcessArguments(value, 5, args))
  return; 

 G4double theta(GetUnit(args[4], "Angle"));
 
 if (theta<=0.)
  return;
 
 G4double phi(theta*G4UIcommand::ConvertToDouble(args[2]));
 G4double delta(theta*G4UIcommand::ConvertToDouble(args[3]));
 theta*=G4UIcommand::ConvertToDouble(args[1]);

 G4ThreeVector axis;
 axis.setRThetaPhi(1., theta/rad, phi/rad);
 detectorLayout->SetPlacementRotation(args[0], G4RotationMatrix(axis, delta/rad));
}



void                                            RadmonDetectorMessenger :: OnSetRelativePlacementPosition(const G4String & value)
{
 G4String args[6];

 if (!ProcessArguments(value, 6, args))
  return; 

 G4double x(GetUnit(args[5], "Length"));
 
 if (x<=0.)
  return;
 
 G4double y(x*G4UIcommand::ConvertToDouble(args[3]));
 G4double z(x*G4UIcommand::ConvertToDouble(args[4]));
 x*=G4UIcommand::ConvertToDouble(args[2]);

 detectorLayout->SetPlacementPosition(args[0], args[1], G4ThreeVector(x, y, z));
}



void                                            RadmonDetectorMessenger :: OnSetRelativePlacementRotation(const G4String & value)
{
 G4String args[6];

 if (!ProcessArguments(value, 6, args))
  return; 

 G4double theta(GetUnit(args[5], "Angle"));
 
 if (theta<=0.)
  return;
 
 G4double phi(theta*G4UIcommand::ConvertToDouble(args[3]));
 G4double delta(theta*G4UIcommand::ConvertToDouble(args[4]));
 theta*=G4UIcommand::ConvertToDouble(args[2]);

 G4ThreeVector axis;
 axis.setRThetaPhi(1., theta, phi);
 detectorLayout->SetPlacementRotation(args[0], args[1], G4RotationMatrix(axis, delta));
}



void                                            RadmonDetectorMessenger :: OnDumpLayout(const G4String & /* value */)
{
 detectorLayout->DumpLayout(G4cout);
 G4cout << G4endl;
}



void                                            RadmonDetectorMessenger :: OnLoad(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 
  
 std::istream * in(OpenForInput(args[0]));
 
 if (!in)
  return;

 if (!detectorLayout->Load(*in)) 
  G4cout << "RadmonDetectorMessenger::OnLoad(): Error reading from file \"" << args[0] << "\"." << G4endl;
  
 delete in;
}



void                                            RadmonDetectorMessenger :: OnSave(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 
  
 std::ostream * out(OpenForOutput(args[0]));
 
 if (!out)
  return;

 if (!detectorLayout->Save(*out))
  G4cout << "RadmonDetectorMessenger::OnSave(): Cannot write layout into file \"" << args[0] << "\"." << G4endl;
  
 delete out;
}
