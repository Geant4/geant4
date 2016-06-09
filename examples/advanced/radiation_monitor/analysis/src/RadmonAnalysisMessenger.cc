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
// File name:     RadmonAnalysisMessenger.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonAnalysisMessenger.cc,v 1.3 2006/06/29 16:07:39 gunter Exp $
// Tag:           $Name: geant4-08-02 $
//

// Messenger commands path
#define COMMANDS_PATH "/radmon/analysis/"

// Include files
#include "RadmonAnalysisMessenger.hh"
#include "RadmonVAnalysisLayout.hh"
#include "G4UnitsTable.hh"



                                                RadmonAnalysisMessenger :: RadmonAnalysisMessenger(RadmonVAnalysisLayout * layout)
:
 RadmonMessenger(COMMANDS_PATH, "Interactive analysis configuration commands."),
 analysisLayout(layout),
 RADMON_INITIALIZE_COMMAND(SetOutputFileName),
 RADMON_INITIALIZE_COMMAND(SetOutputFileFormat),
 RADMON_INITIALIZE_COMMAND(CreateSensitiveDetector),
 RADMON_INITIALIZE_COMMAND(SetSensitiveDetectorType),
 RADMON_INITIALIZE_COMMAND(RemoveSensitiveDetector),
 RADMON_INITIALIZE_COMMAND(CreateSensitiveDetectorType),
 RADMON_INITIALIZE_COMMAND(RemoveSensitiveDetectorType),
 RADMON_INITIALIZE_COMMAND(AppendDataAnalysisToSensitiveDetectorType),
 RADMON_INITIALIZE_COMMAND(SetDataAnalysisType),
 RADMON_INITIALIZE_COMMAND(RemoveDataAnalysis),
 RADMON_INITIALIZE_COMMAND(SetDataAnalysisAttribute),
 RADMON_INITIALIZE_COMMAND(ClearDataAnalysisAttribute),
 RADMON_INITIALIZE_COMMAND(DumpLayout),
 RADMON_INITIALIZE_COMMAND(Load),
 RADMON_INITIALIZE_COMMAND(Save)
{
 if (layout==0)
  G4Exception("RadmonAnalysisMessenger::RadmonAnalysisMessenger: layout==0.");

 RADMON_CREATE_COMMAND_1ARG (SetOutputFileName,                         "Sets the name of the analysis file",                         "fileName");
 RADMON_CREATE_COMMAND_1ARG (SetOutputFileFormat,                       "Sets the format of the analysis file",                       "format");
 RADMON_CREATE_COMMAND_2ARGS(CreateSensitiveDetector,                   "Creates a sensitive detector",                               "sensitiveDetectorName", "sensitiveDetectorType");
 RADMON_CREATE_COMMAND_2ARGS(SetSensitiveDetectorType,                  "Changes the type of sensitive detector",                     "sensitiveDetectorName", "sensitiveDetectorType");
 RADMON_CREATE_COMMAND_1ARG (RemoveSensitiveDetector,                   "Removes a sensitive detector",                               "sensitiveDetectorName");
 RADMON_CREATE_COMMAND_1ARG (CreateSensitiveDetectorType,               "Creates a sensitive detector type",                          "sensitiveDetectorType");
 RADMON_CREATE_COMMAND_1ARG (RemoveSensitiveDetectorType,               "Removes a sensitive detector type",                          "sensitiveDetectorType");
 RADMON_CREATE_COMMAND_2ARGS(AppendDataAnalysisToSensitiveDetectorType, "Appends a analysis module to the sensitive detector type",   "sensitiveDetectorType", "analysisModuleLabel");
 RADMON_CREATE_COMMAND_3ARGS(SetDataAnalysisType,                       "Sets the type of analysis module",                           "sensitiveDetectorType", "analysisModuleLabel", "analysisModuleType");
 RADMON_CREATE_COMMAND_2ARGS(RemoveDataAnalysis,                        "Removes a analysis module from the sensitive detector type", "sensitiveDetectorType", "analysisModuleLabel");
 RADMON_CREATE_COMMAND_4ARGS(SetDataAnalysisAttribute,                  "Sets an attribute of a analysis module",                     "sensitiveDetectorType", "analysisModuleLabel", "attributeName",     "attributeValue");
 RADMON_CREATE_COMMAND_3ARGS(ClearDataAnalysisAttribute,                "Clears an attribute of a analysis module",                   "sensitiveDetectorType", "analysisModuleLabel", "attributeName");
 RADMON_CREATE_COMMAND_0ARGS(DumpLayout,                                "Dumps the layout");
 RADMON_CREATE_COMMAND_1ARG (Load,                                      "Loads the configuration from file",                          "fileName");
 RADMON_CREATE_COMMAND_1ARG (Save,                                      "Saves the configuration from file",                          "fileName");
}



                                                RadmonAnalysisMessenger :: ~RadmonAnalysisMessenger()
{
 RADMON_DESTROY_COMMAND(Save);
 RADMON_DESTROY_COMMAND(Load);
 RADMON_DESTROY_COMMAND(DumpLayout);
 RADMON_DESTROY_COMMAND(ClearDataAnalysisAttribute);
 RADMON_DESTROY_COMMAND(SetDataAnalysisAttribute);
 RADMON_DESTROY_COMMAND(RemoveDataAnalysis);
 RADMON_DESTROY_COMMAND(SetDataAnalysisType);
 RADMON_DESTROY_COMMAND(AppendDataAnalysisToSensitiveDetectorType);
 RADMON_DESTROY_COMMAND(RemoveSensitiveDetectorType);
 RADMON_DESTROY_COMMAND(CreateSensitiveDetectorType);
 RADMON_DESTROY_COMMAND(RemoveSensitiveDetector);
 RADMON_DESTROY_COMMAND(SetSensitiveDetectorType);
 RADMON_DESTROY_COMMAND(CreateSensitiveDetector);
 RADMON_DESTROY_COMMAND(SetOutputFileFormat);
 RADMON_DESTROY_COMMAND(SetOutputFileName);
}





G4String                                        RadmonAnalysisMessenger :: GetCurrentValue(G4UIcommand * /* command */)
{
 G4cout << "RadmonAnalysisMessenger::GetCurrentValue(): Not supported" << G4endl;
 
 return G4String();
}



void                                            RadmonAnalysisMessenger :: SetNewValue(G4UIcommand * command, G4String newValue)
{
 RADMON_BEGIN_LIST_SET_COMMANDS
  RADMON_SET_COMMAND(SetOutputFileName)
  RADMON_SET_COMMAND(SetOutputFileFormat)
  RADMON_SET_COMMAND(CreateSensitiveDetector)
  RADMON_SET_COMMAND(SetSensitiveDetectorType)
  RADMON_SET_COMMAND(RemoveSensitiveDetector)
  RADMON_SET_COMMAND(CreateSensitiveDetectorType)
  RADMON_SET_COMMAND(RemoveSensitiveDetectorType)
  RADMON_SET_COMMAND(AppendDataAnalysisToSensitiveDetectorType)
  RADMON_SET_COMMAND(SetDataAnalysisType)
  RADMON_SET_COMMAND(RemoveDataAnalysis)
  RADMON_SET_COMMAND(SetDataAnalysisAttribute)
  RADMON_SET_COMMAND(ClearDataAnalysisAttribute)
  RADMON_SET_COMMAND(DumpLayout)
  RADMON_SET_COMMAND(Load)
  RADMON_SET_COMMAND(Save)
 RADMON_END_LIST_SET_COMMANDS
}





// Events
void                                            RadmonAnalysisMessenger :: OnSetOutputFileName(const G4String & value)
{
 G4String args;

 if (!ProcessArguments(value, 1, &args))
  return; 
  
 analysisLayout->SetOutputFileName(args);
}



void                                            RadmonAnalysisMessenger :: OnSetOutputFileFormat(const G4String & value)
{
 G4String args;

 if (!ProcessArguments(value, 1, &args))
  return; 
  
 analysisLayout->SetOutputFileFormat(args);
}



void                                            RadmonAnalysisMessenger :: OnCreateSensitiveDetector(const G4String & value)
{
 G4String args[2];

 if (!ProcessArguments(value, 2, args))
  return; 
  
 analysisLayout->CreateSensitiveDetector(args[0], args[1]);
}



void                                            RadmonAnalysisMessenger :: OnSetSensitiveDetectorType(const G4String & value)
{
 G4String args[2];

 if (!ProcessArguments(value, 2, args))
  return; 
  
 analysisLayout->SetSensitiveDetectorType(args[0], args[1]);
}



void                                            RadmonAnalysisMessenger :: OnRemoveSensitiveDetector(const G4String & value)
{
 G4String args;

 if (!ProcessArguments(value, 1, &args))
  return; 
  
 analysisLayout->RemoveSensitiveDetector(args);
}



void                                            RadmonAnalysisMessenger :: OnCreateSensitiveDetectorType(const G4String & value)
{
 G4String args;

 if (!ProcessArguments(value, 1, &args))
  return; 
  
 analysisLayout->CreateSensitiveDetectorType(args);
}



void                                            RadmonAnalysisMessenger :: OnRemoveSensitiveDetectorType(const G4String & value)
{
 G4String args;

 if (!ProcessArguments(value, 1, &args))
  return; 
  
 analysisLayout->RemoveSensitiveDetectorType(args);
}



void                                            RadmonAnalysisMessenger :: OnAppendDataAnalysisToSensitiveDetectorType(const G4String & value)
{
 G4String args[2];

 if (!ProcessArguments(value, 2, args))
  return; 
  
 analysisLayout->AppendDataAnalysisToSensitiveDetectorType(args[0], args[1]);
}



void                                            RadmonAnalysisMessenger :: OnSetDataAnalysisType(const G4String & value)
{
 G4String args[3];

 if (!ProcessArguments(value, 3, args))
  return; 
  
 analysisLayout->SetDataAnalysisType(args[0], args[1], args[2]);
}



void                                            RadmonAnalysisMessenger :: OnRemoveDataAnalysis(const G4String & value)
{
 G4String args[2];

 if (!ProcessArguments(value, 2, args))
  return; 
  
 analysisLayout->RemoveDataAnalysis(args[0], args[1]);
}



void                                            RadmonAnalysisMessenger :: OnSetDataAnalysisAttribute(const G4String & value)
{
 G4String args[4];

 if (!ProcessArguments(value, 4, args))
  return; 
  
 analysisLayout->SetDataAnalysisAttribute(args[0], args[1], args[2], args[3]);
}



void                                            RadmonAnalysisMessenger :: OnClearDataAnalysisAttribute(const G4String & value)
{
 G4String args[3];

 if (!ProcessArguments(value, 3, args))
  return; 
  
 analysisLayout->ClearDataAnalysisAttribute(args[0], args[1], args[2]);
}



void                                            RadmonAnalysisMessenger :: OnDumpLayout(const G4String & /* value */)
{
 analysisLayout->DumpLayout(G4cout);
 G4cout << G4endl;
}



void                                            RadmonAnalysisMessenger :: OnLoad(const G4String & value)
{
 G4String args;

 if (!ProcessArguments(value, 1, &args))
  return; 
  
 std::istream * in(OpenForInput(args));
 
 if (!in)
  return;

 if (!analysisLayout->Load(*in)) 
  G4cout << "RadmonAnalysisMessenger::OnLoad(): Error reading from file \"" << args << "\"." << G4endl;
  
 delete in;
}



void                                            RadmonAnalysisMessenger :: OnSave(const G4String & value)
{
 G4String args;

 if (!ProcessArguments(value, 1, &args))
  return; 
  
 std::ostream * out(OpenForOutput(args));
 
 if (!out)
  return;

 if (!analysisLayout->Save(*out))
  G4cout << "RadmonAnalysisMessenger::OnSave(): Cannot write layout into file \"" << args << "\"." << G4endl;
  
 delete out;
}
