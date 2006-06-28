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
// File name:     RadmonPhysicsMessenger.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsMessenger.cc,v 1.3 2006-06-28 13:56:25 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Messenger commands path
#define COMMANDS_PATH "/radmon/physics/"

// Include files
#include "RadmonPhysicsMessenger.hh"
#include "RadmonVPhysicsLayout.hh"
#include "G4UnitsTable.hh"



                                                RadmonPhysicsMessenger :: RadmonPhysicsMessenger(RadmonVPhysicsLayout * layout)
:
 RadmonMessenger(COMMANDS_PATH, "Interactive physics list commands."),
 physicsLayout(layout),
 RADMON_INITIALIZE_COMMAND(AddPhysicsList),
 RADMON_INITIALIZE_COMMAND(RemovePhysicsList),
 RADMON_INITIALIZE_COMMAND(SetPhysicsListAttribute),
 RADMON_INITIALIZE_COMMAND(ClearPhysicsListAttribute),
 RADMON_INITIALIZE_COMMAND(DumpLayout),
 RADMON_INITIALIZE_COMMAND(Load),
 RADMON_INITIALIZE_COMMAND(Save)
{
 if (layout==0)
  G4Exception("RadmonPhysicsMessenger::RadmonPhysicsMessenger: layout==0.");

 RADMON_CREATE_COMMAND_1ARG (AddPhysicsList,                    "Enables a group of physics processes",    "physiscsListName");
 RADMON_CREATE_COMMAND_1ARG (RemovePhysicsList,                 "Disables a group of physics processes",   "physicsListName");
 RADMON_CREATE_COMMAND_3ARGS(SetPhysicsListAttribute,           "Set the attributes of a physics list",    "physicsListName", "attributeName", "value");
 RADMON_CREATE_COMMAND_2ARGS(ClearPhysicsListAttribute,         "Remove the attributes of a physics list", "physicsListName", "attributeName");
 RADMON_CREATE_COMMAND_0ARGS(DumpLayout,                        "Print out the current layout");
 RADMON_CREATE_COMMAND_1ARG (Load,                              "Loads a layout from file",                "fileName");
 RADMON_CREATE_COMMAND_1ARG (Save,                              "Saves a layout to file",                  "fileName");
}



                                                RadmonPhysicsMessenger :: ~RadmonPhysicsMessenger()
{
 RADMON_DESTROY_COMMAND(Save);
 RADMON_DESTROY_COMMAND(Load);
 RADMON_DESTROY_COMMAND(DumpLayout);
 RADMON_DESTROY_COMMAND(ClearPhysicsListAttribute);
 RADMON_DESTROY_COMMAND(SetPhysicsListAttribute);
 RADMON_DESTROY_COMMAND(RemovePhysicsList);
 RADMON_DESTROY_COMMAND(AddPhysicsList);
}





G4String                                        RadmonPhysicsMessenger :: GetCurrentValue(G4UIcommand * /* command */)
{
 G4cout << "RadmonPhysicsMessenger::GetCurrentValue(): Not supported" << G4endl;
 
 return G4String();
}



void                                            RadmonPhysicsMessenger :: SetNewValue(G4UIcommand * command, G4String newValue)
{
 RADMON_BEGIN_LIST_SET_COMMANDS
  RADMON_SET_COMMAND(AddPhysicsList)
  RADMON_SET_COMMAND(RemovePhysicsList)
  RADMON_SET_COMMAND(SetPhysicsListAttribute)
  RADMON_SET_COMMAND(ClearPhysicsListAttribute)
  RADMON_SET_COMMAND(DumpLayout)
  RADMON_SET_COMMAND(Load)
  RADMON_SET_COMMAND(Save)
 RADMON_END_LIST_SET_COMMANDS
}





// Events
void                                            RadmonPhysicsMessenger :: OnAddPhysicsList(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 
 
 physicsLayout->AddPhysicsList(args[0]); 
}



void                                            RadmonPhysicsMessenger :: OnRemovePhysicsList(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 
 
 physicsLayout->RemovePhysicsList(args[0]); 
}



void                                            RadmonPhysicsMessenger :: OnSetPhysicsListAttribute(const G4String & value)
{
 G4String args[3];

 if (!ProcessArguments(value, 3, args))
  return; 
 
 physicsLayout->SetPhysicsListAttribute(args[0], args[1], args[2]); 
}



void                                            RadmonPhysicsMessenger :: OnClearPhysicsListAttribute(const G4String & value)
{
 G4String args[2];

 if (!ProcessArguments(value, 2, args))
  return; 
 
 physicsLayout->ClearPhysicsListAttribute(args[0], args[1]);
}



void                                            RadmonPhysicsMessenger :: OnDumpLayout(const G4String & /* value */)
{
 physicsLayout->DumpLayout(G4cout);
 G4cout << G4endl;
}



void                                            RadmonPhysicsMessenger :: OnLoad(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 
  
 std::istream * in(OpenForInput(args[0]));
 
 if (!in)
  return;

 if (!physicsLayout->Load(*in)) 
  G4cout << "RadmonPhysicsMessenger::OnLoad(): Error reading from file \"" << args[0] << "\"." << G4endl;
  
 delete in;
}



void                                            RadmonPhysicsMessenger :: OnSave(const G4String & value)
{
 G4String args[1];

 if (!ProcessArguments(value, 1, args))
  return; 
  
 std::ostream * out(OpenForOutput(args[0]));
 
 if (!out)
  return;

 if (!physicsLayout->Save(*out))
  G4cout << "RadmonPhysicsMessenger::OnSave(): Cannot write layout into file \"" << args[0] << "\"." << G4endl;
  
 delete out;
}
