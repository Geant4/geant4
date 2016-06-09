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
// File name:     RadmonGeneratorMessenger.cc
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorMessenger.cc,v 1.3 2006/06/29 16:16:27 gunter Exp $
// Tag:           $Name: geant4-09-00 $
//

// Messenger commands path
#define COMMANDS_PATH "/radmon/generator/"

// Include files
#include "RadmonGeneratorMessenger.hh"
#include "RadmonVGeneratorLayout.hh"



                                                RadmonGeneratorMessenger :: RadmonGeneratorMessenger(RadmonVGeneratorLayout * layout)
:
 RadmonMessenger(COMMANDS_PATH, "Interactive beam generation commands."),
 generatorLayout(layout),
 RADMON_INITIALIZE_COMMAND(InsertSource),
 RADMON_INITIALIZE_COMMAND(SetRelativeSourceIntensity),
 RADMON_INITIALIZE_COMMAND(RemoveSource),
 RADMON_INITIALIZE_COMMAND(AppendSourceAlgorithm),
 RADMON_INITIALIZE_COMMAND(SetSourceAlgorithmType),
 RADMON_INITIALIZE_COMMAND(RemoveSourceAlgorithm),
 RADMON_INITIALIZE_COMMAND(SetSourceAlgorithmAttribute),
 RADMON_INITIALIZE_COMMAND(ClearSourceAlgorithmAttribute),
 RADMON_INITIALIZE_COMMAND(Load),
 RADMON_INITIALIZE_COMMAND(Save),
 RADMON_INITIALIZE_COMMAND(DumpLayout)
{
 if (layout==0)
  G4Exception("RadmonGeneratorMessenger::RadmonGeneratorMessenger: layout==0.");

 RADMON_CREATE_COMMAND_1ARG (InsertSource,                  "Inserts a new source",                                "sourcename");
 RADMON_CREATE_COMMAND_2ARGS(SetRelativeSourceIntensity,    "Sets the relative intensity of the source",           "sourcename", "intensity");
 RADMON_CREATE_COMMAND_1ARG (RemoveSource,                  "Removes a source from the layout",                    "sourcename");
 RADMON_CREATE_COMMAND_2ARGS(AppendSourceAlgorithm,         "Appends an algorithm to a source",                    "sourcename", "algorithm");
 RADMON_CREATE_COMMAND_3ARGS(SetSourceAlgorithmType,        "Sets the algorithm type of an algorithm or a source", "sourcename", "algorithm", "typename");
 RADMON_CREATE_COMMAND_2ARGS(RemoveSourceAlgorithm,         "Removes an algorithm from a source",                  "sourcename", "algorithm");
 RADMON_CREATE_COMMAND_4ARGS(SetSourceAlgorithmAttribute,   "Sets an attribute of the algorithm of a source",      "sourcename", "algorithm", "attributename", "value");
 RADMON_CREATE_COMMAND_3ARGS(ClearSourceAlgorithmAttribute, "Removes an attribute from the algorithm of a source", "sourcename", "algorithm", "attributename");
 RADMON_CREATE_COMMAND_1ARG (Load,                          "Loads a layout from file",                            "filename");
 RADMON_CREATE_COMMAND_1ARG (Save,                          "Saves a layout to file",                              "filename");
 RADMON_CREATE_COMMAND_0ARGS(DumpLayout,                    "Print out the current layout");
}



                                                RadmonGeneratorMessenger :: ~RadmonGeneratorMessenger()
{
 RADMON_DESTROY_COMMAND(DumpLayout);
 RADMON_DESTROY_COMMAND(Save);
 RADMON_DESTROY_COMMAND(Load);
 RADMON_DESTROY_COMMAND(ClearSourceAlgorithmAttribute);
 RADMON_DESTROY_COMMAND(SetSourceAlgorithmAttribute);
 RADMON_DESTROY_COMMAND(RemoveSourceAlgorithm);
 RADMON_DESTROY_COMMAND(SetSourceAlgorithmType);
 RADMON_DESTROY_COMMAND(AppendSourceAlgorithm);
 RADMON_DESTROY_COMMAND(RemoveSource);
 RADMON_DESTROY_COMMAND(SetRelativeSourceIntensity);
 RADMON_DESTROY_COMMAND(InsertSource);
}





G4String                                        RadmonGeneratorMessenger :: GetCurrentValue(G4UIcommand * /* command */)
{
 G4cout << "RadmonGeneratorMessenger::GetCurrentValue(): Not supported" << G4endl;
 
 return G4String();
}



void                                            RadmonGeneratorMessenger :: SetNewValue(G4UIcommand * command, G4String newValue)
{
 RADMON_BEGIN_LIST_SET_COMMANDS
  RADMON_SET_COMMAND(InsertSource)
  RADMON_SET_COMMAND(SetRelativeSourceIntensity)
  RADMON_SET_COMMAND(RemoveSource)
  RADMON_SET_COMMAND(AppendSourceAlgorithm)
  RADMON_SET_COMMAND(SetSourceAlgorithmType)
  RADMON_SET_COMMAND(RemoveSourceAlgorithm)
  RADMON_SET_COMMAND(SetSourceAlgorithmAttribute)
  RADMON_SET_COMMAND(ClearSourceAlgorithmAttribute)
  RADMON_SET_COMMAND(Load)
  RADMON_SET_COMMAND(Save)
  RADMON_SET_COMMAND(DumpLayout)
 RADMON_END_LIST_SET_COMMANDS
}





// Events
void                                            RadmonGeneratorMessenger :: OnInsertSource(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 
 
 generatorLayout->InsertSource(args[0]);
}



void                                            RadmonGeneratorMessenger :: OnSetRelativeSourceIntensity(const G4String & value)
{
 G4String args[2];

 if (!ProcessArguments(value, 2, args))
  return; 
 
 G4double intensity(G4UIcommand::ConvertToDouble(args[1]));
 
 if (intensity<=0.)
 {
  G4cout << "RadmonGeneratorMessenger::OnSetRelativeSourceIntensity(): insensitymust be greater than 0." << G4endl;
  return;  
 }
 
 generatorLayout->SetRelativeSourceIntensity(args[0], intensity);
}



void                                            RadmonGeneratorMessenger :: OnRemoveSource(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 
 
 generatorLayout->RemoveSource(args[0]);
}





void                                            RadmonGeneratorMessenger :: OnAppendSourceAlgorithm(const G4String & value)
{
 G4String args[2];

 if (!ProcessArguments(value, 2, args))
  return; 
 
 generatorLayout->AppendSourceAlgorithm(args[0], args[1]);
}



void                                            RadmonGeneratorMessenger :: OnSetSourceAlgorithmType(const G4String & value)
{
 G4String args[3];

 if (!ProcessArguments(value, 3, args))
  return; 
 
 generatorLayout->SetSourceAlgorithmType(args[0], args[1], args[2]);
}



void                                            RadmonGeneratorMessenger :: OnRemoveSourceAlgorithm(const G4String & value)
{
 G4String args[2];

 if (!ProcessArguments(value, 2, args))
  return; 
 
 generatorLayout->RemoveSourceAlgorithm(args[0], args[1]);
}





void                                            RadmonGeneratorMessenger :: OnSetSourceAlgorithmAttribute(const G4String & value)
{
 G4String args[4];

 if (!ProcessArguments(value, 4, args))
  return; 
 
 generatorLayout->SetSourceAlgorithmAttribute(args[0], args[1], args[2], args[3]);
}



void                                            RadmonGeneratorMessenger :: OnClearSourceAlgorithmAttribute(const G4String & value)
{
 G4String args[3];

 if (!ProcessArguments(value, 3, args))
  return; 
 
 generatorLayout->ClearSourceAlgorithmAttribute(args[0], args[1], args[2]);
}





void                                            RadmonGeneratorMessenger :: OnLoad(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 
  
 std::istream * in(OpenForInput(args[0]));
 
 if (!in)
  return;

 if (!generatorLayout->Load(*in)) 
  G4cout << "RadmonGeneratorMessenger::OnLoad(): Error reading from file \"" << args[0] << "\"." << G4endl;
  
 delete in;
}



void                                            RadmonGeneratorMessenger :: OnSave(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 
  
 std::ostream * out(OpenForOutput(args[0]));
 
 if (!out)
  return;

 if (!generatorLayout->Save(*out))
  G4cout << "RadmonGeneratorMessenger::OnSave(): Cannot write layout into file \"" << args[0] << "\"." << G4endl;
  
 delete out;
}



void                                            RadmonGeneratorMessenger :: OnDumpLayout(const G4String & /* value */)
{
 generatorLayout->DumpLayout(G4cout);
 G4cout << G4endl;
}
