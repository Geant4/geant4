//
// File name:     RadmonDetectorMessenger.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMessenger.cc,v 1.2 2005-09-14 12:28:31 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorMessenger.hh"
#include "RadmonVDetectorLayout.hh"

#include "globals.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "RadmonTokenizer.hh"

#include <fstream>

#define COMMANDS_PATH "/radmon/detector/"

#define INITIALIZE_COMMAND(command)             cmd ## command (0)

#define CREATE_COMMAND(command, guidance)       cmd ## command = new G4UIcommand(COMMANDS_PATH #command , this);                                                 \
                                                                                                                                                                 \
                                                if (! cmd ## command )                                                                                           \
                                                {                                                                                                                \
                                                 G4cerr << "RadmonDetectorMessenger::RadmonDetectorMessenger: Command \"" #command "\" not allocated."<< G4endl; \
                                                 return;                                                                                                         \
                                                }                                                                                                                \
                                                                                                                                                                 \
                                                cmd ## command ->SetGuidance(guidance)

#define DESTROY_COMMAND(command)                if ( cmd ## command )   \
                                                 delete cmd ## command;
                                                
#define _ADD_ARG(command, argName)              cmd ## command ->SetParameter(new G4UIparameter(argName, 's', false))

#define CREATE_COMMAND_1ARG(command, guidance, argName0)                           \
                                                CREATE_COMMAND(command, guidance); \
                                                _ADD_ARG(command, argName0)

#define CREATE_COMMAND_2ARG(command, guidance, argName0, argName1)                                \
                                                CREATE_COMMAND_1ARG(command, guidance, argName0); \
                                                _ADD_ARG(command, argName1)

#define CREATE_COMMAND_3ARG(command, guidance, argName0, argName1, argName2)                                \
                                                CREATE_COMMAND_2ARG(command, guidance, argName0, argName1); \
                                                _ADD_ARG(command, argName2)

#define CREATE_COMMAND_4ARG(command, guidance, argName0, argName1, argName2, argName3)                                \
                                                CREATE_COMMAND_3ARG(command, guidance, argName0, argName1, argName2); \
                                                _ADD_ARG(command, argName3)

#define CREATE_COMMAND_5ARG(command, guidance, argName0, argName1, argName2, argName3, argName4)                                \
                                                CREATE_COMMAND_4ARG(command, guidance, argName0, argName1, argName2, argName3); \
                                                _ADD_ARG(command, argName4)

#define CREATE_COMMAND_6ARG(command, guidance, argName0, argName1, argName2, argName3, argName4, argName5)                                \
                                                CREATE_COMMAND_5ARG(command, guidance, argName0, argName1, argName2, argName3, argName4); \
                                                _ADD_ARG(command, argName5)



                                                RadmonDetectorMessenger :: RadmonDetectorMessenger(RadmonVDetectorLayout * layout)
:
 detectorLayout(layout),
 directory(0),
 INITIALIZE_COMMAND(EnableEnvironment),
 INITIALIZE_COMMAND(DisableEnvironment),
 INITIALIZE_COMMAND(SetEnvironmentType),
 INITIALIZE_COMMAND(SetEnvironmentAttribute),
 INITIALIZE_COMMAND(ClearEnvironmentAttribute),
 INITIALIZE_COMMAND(CreateMultilayer),
 INITIALIZE_COMMAND(RemoveMultilayer),
 INITIALIZE_COMMAND(SetMultilayerWidth),
 INITIALIZE_COMMAND(SetMultilayerHeight),
 INITIALIZE_COMMAND(AppendLayerToMultilayer),
 INITIALIZE_COMMAND(RemoveLayerFromMultilayer),
 INITIALIZE_COMMAND(RemoveAllLayersFromMultilayer),
 INITIALIZE_COMMAND(SetLayerThickness),
 INITIALIZE_COMMAND(SetLayerType),
 INITIALIZE_COMMAND(SetLayerAttribute),
 INITIALIZE_COMMAND(ClearLayerAttribute),
 INITIALIZE_COMMAND(CreatePlacement),
 INITIALIZE_COMMAND(RemovePlacement),
 INITIALIZE_COMMAND(SetPlacementPosition),
 INITIALIZE_COMMAND(SetPlacementRotation),
 INITIALIZE_COMMAND(SetRelativePlacementPosition),
 INITIALIZE_COMMAND(SetRelativePlacementRotation),
 INITIALIZE_COMMAND(DumpLayout),
 INITIALIZE_COMMAND(Load),
 INITIALIZE_COMMAND(Save)
{
 if (layout==0)
  G4Exception("RadmonDetectorMessenger::RadmonDetectorMessenger: layout==0.");
 
 directory=new G4UIdirectory(COMMANDS_PATH);
 
 if (directory==0)
 {
  G4cerr << "RadmonDetectorMessenger::RadmonDetectorMessenger: \"" COMMANDS_PATH "\" directory not allocated." << G4endl;
  return;
 }
 
 directory->SetGuidance("Intercative detector construction commands.");

 CREATE_COMMAND(EnableEnvironment, "Enables the environment");
 CREATE_COMMAND(DisableEnvironment, "Disables the environment");
 CREATE_COMMAND_1ARG(SetEnvironmentType, "Define the type of environment", "type");
 CREATE_COMMAND_2ARG(SetEnvironmentAttribute, "Define an attribute of the environment", "name", "value");
 CREATE_COMMAND_1ARG(ClearEnvironmentAttribute, "Removes an attribute of the environment", "name"); 
 CREATE_COMMAND_1ARG(CreateMultilayer, "Creates a new multilayer", "name");
 CREATE_COMMAND_1ARG(RemoveMultilayer, "Removes a multilayer", "name");
 CREATE_COMMAND_3ARG(SetMultilayerWidth, "Set the width of a multilayer", "name", "width", "unit");
 CREATE_COMMAND_3ARG(SetMultilayerHeight, "Set the height of a multilayer", "name", "height", "unit");
 CREATE_COMMAND_2ARG(AppendLayerToMultilayer, "Adds a layer to a multilayer", "multilayerName", "layerName");
 CREATE_COMMAND_2ARG(RemoveLayerFromMultilayer, "Removes a layer from a multilayer", "multilayerName", "layerName");
 CREATE_COMMAND_1ARG(RemoveAllLayersFromMultilayer, "Removes all the layers of a multilayer", "multilayerName");
 CREATE_COMMAND_4ARG(SetLayerThickness, "Set the thickness of the layer of a multilayer", "multilayerName", "layerName", "width", "unit");
 CREATE_COMMAND_3ARG(SetLayerType, "Set the layer type of a multilayer", "multilayerName", "layerName", "type");
 CREATE_COMMAND_4ARG(SetLayerAttribute, "Set an attribute of a layer of a multilayer", "multilayerName", "layerName", "attributeName", "value");
 CREATE_COMMAND_3ARG(ClearLayerAttribute, "Removes an attribute of a layer of a multilayer", "multilayerName", "layerName", "attributeName");
 CREATE_COMMAND_2ARG(CreatePlacement, "Creates a new multilayer placement", "placementName", "multilayerName");
 CREATE_COMMAND_1ARG(RemovePlacement, "Removes a multilayer placement", "name");
 CREATE_COMMAND_5ARG(SetPlacementPosition,  "Changes the position of a placement", "name", "x", "y", "z", "unit");
 CREATE_COMMAND_5ARG(SetPlacementRotation,  "Changes the orientation of a placement performing a rotation of delta along axis (theta, phi)", "name", "theta", "phi", "delta", "unit");
 CREATE_COMMAND_6ARG(SetRelativePlacementPosition,  "Changes the position of a placement relative to another placement", "name", "fromName", "x", "y", "z", "unit");
 CREATE_COMMAND_6ARG(SetRelativePlacementRotation,  "Changes the orientation of a placement relative to another placement", "name", "fromName", "theta", "phi", "delta", "unit");
 CREATE_COMMAND(DumpLayout, "Print out the current layout");
 CREATE_COMMAND_1ARG(Load, "Loads a layout from file", "fileName");
 CREATE_COMMAND_1ARG(Save, "Saves a layout to file", "fileName");
}



                                                RadmonDetectorMessenger :: ~RadmonDetectorMessenger()
{
 DESTROY_COMMAND(Save);
 DESTROY_COMMAND(Load);
 DESTROY_COMMAND(DumpLayout);
 DESTROY_COMMAND(SetRelativePlacementRotation);
 DESTROY_COMMAND(SetRelativePlacementPosition);
 DESTROY_COMMAND(SetPlacementRotation);
 DESTROY_COMMAND(SetPlacementPosition);
 DESTROY_COMMAND(RemovePlacement);
 DESTROY_COMMAND(CreatePlacement);
 DESTROY_COMMAND(ClearLayerAttribute);
 DESTROY_COMMAND(SetLayerAttribute);
 DESTROY_COMMAND(SetLayerType);
 DESTROY_COMMAND(SetLayerThickness);
 DESTROY_COMMAND(RemoveAllLayersFromMultilayer);
 DESTROY_COMMAND(RemoveLayerFromMultilayer);
 DESTROY_COMMAND(AppendLayerToMultilayer);
 DESTROY_COMMAND(SetMultilayerHeight);
 DESTROY_COMMAND(SetMultilayerWidth);
 DESTROY_COMMAND(RemoveMultilayer);
 DESTROY_COMMAND(CreateMultilayer);
 DESTROY_COMMAND(ClearEnvironmentAttribute);
 DESTROY_COMMAND(SetEnvironmentAttribute);
 DESTROY_COMMAND(SetEnvironmentType);
 DESTROY_COMMAND(DisableEnvironment);
 DESTROY_COMMAND(EnableEnvironment);

 if (directory)
  delete directory;
}





G4String                                        RadmonDetectorMessenger :: GetCurrentValue(G4UIcommand * /* command */)
{
 G4cerr << "RadmonDetectorMessenger::GetCurrentValue(): Not supported" << G4endl;
 
 return G4String();
}



#define BEGIN_LIST_SET_COMMANDS                 if (!command)                                                             \
                                                 G4cerr << "RadmonDetectorMessenger::SetNewValue: command==0." << G4endl; \
                                                else

#define  SET_COMMAND(name)                      if ( cmd ## name == command) \
                                                 On ## name (newValue);      \
                                                else

#define END_LIST_SET_COMMANDS                   G4cerr << "RadmonDetectorMessenger::SetNewValue: Command \"" << command->GetCommandPath() << "\" not supported." << G4endl;

void                                            RadmonDetectorMessenger :: SetNewValue(G4UIcommand * command, G4String newValue)
{
 BEGIN_LIST_SET_COMMANDS
  SET_COMMAND(EnableEnvironment)
  SET_COMMAND(DisableEnvironment)
  SET_COMMAND(SetEnvironmentType)
  SET_COMMAND(SetEnvironmentAttribute)
  SET_COMMAND(ClearEnvironmentAttribute)
  SET_COMMAND(CreateMultilayer)
  SET_COMMAND(RemoveMultilayer)
  SET_COMMAND(SetMultilayerWidth)
  SET_COMMAND(SetMultilayerHeight)
  SET_COMMAND(AppendLayerToMultilayer)
  SET_COMMAND(RemoveLayerFromMultilayer)
  SET_COMMAND(RemoveAllLayersFromMultilayer)
  SET_COMMAND(SetLayerThickness)
  SET_COMMAND(SetLayerType)
  SET_COMMAND(SetLayerAttribute)
  SET_COMMAND(ClearLayerAttribute)
  SET_COMMAND(CreatePlacement)
  SET_COMMAND(RemovePlacement)
  SET_COMMAND(SetPlacementPosition)
  SET_COMMAND(SetPlacementRotation)
  SET_COMMAND(SetRelativePlacementPosition)
  SET_COMMAND(SetRelativePlacementRotation)
  SET_COMMAND(DumpLayout)
  SET_COMMAND(Load)
  SET_COMMAND(Save)
 END_LIST_SET_COMMANDS
}



G4bool                                          RadmonDetectorMessenger :: ProcessArguments(const G4String & rawArguments, G4int nArgs, G4String * arguments) const
{
 RadmonTokenizer args(rawArguments);
 
 for (G4int i(0); i<nArgs; i++)
 {
  if (args.eos())
  {
   G4cout << "RadmonDetectorMessenger::ProcessArguments: " << (nArgs-i) <<  " arguments missing." << G4endl;
   return false;
  }

  arguments[i]=args();
 }

 if (!args.eos())
 {
  G4cout << "RadmonDetectorMessenger::ProcessArguments: Unexpected arguments after \"" << arguments[nArgs-1] << "\"." << G4endl;
  return false;
 }
 
 return true;
}



G4double                                        RadmonDetectorMessenger :: GetUnit(const G4String & unitStr, const char * category) const
{
 G4double dbl(G4UnitDefinition::GetValueOf(unitStr));

 if (dbl<=0.)
 {
  G4cout << "RadmonDetectorMessenger::GetUnit(): Unknown unit of measure \"" << unitStr << "\"." << G4endl;
  return -1.;
 }
 
 if (G4UnitDefinition::GetCategory(unitStr)!=category)
 {
  G4cout << "RadmonDetectorMessenger::GetUnit(): Unit of measure \"" << unitStr << "\" is not a " << category << '.' << G4endl;
  return -1.;
 }
 
 return dbl;
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

 std::ifstream in(args[0].data());

 if (!in.good())
 {
  G4cout << "RadmonDetectorMessenger::OnLoad(): Cannot open file \"" << args[0] << "\" for read." << G4endl;
  return;
 }

 if (!detectorLayout->Load(in)) 
  G4cout << "RadmonDetectorMessenger::OnLoad(): Error reading from file \"" << args[0] << "\"." << G4endl;
}



void                                            RadmonDetectorMessenger :: OnSave(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 

 std::ofstream out(args[0].data());

 if (!out.good())
 {
  G4cout << "RadmonDetectorMessenger::OnSave(): Cannot open file \"" << args[0] << "\" for write." << G4endl;
  return;
 }

 if (!detectorLayout->Save(out))
  G4cout << "RadmonDetectorMessenger::OnSave(): Cannot write layout into file \"" << args[0] << "\"." << G4endl;
}
