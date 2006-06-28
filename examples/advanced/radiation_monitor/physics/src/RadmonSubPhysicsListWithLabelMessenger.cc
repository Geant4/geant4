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
// File name:     RadmonSubPhysicsListWithLabelMessenger.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSubPhysicsListWithLabelMessenger.cc,v 1.4 2006-06-28 13:56:51 gunter Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Messenger commands path
#define COMMANDS_PATH "/radmon/subPhysicsListFactory/"

// Include files
#include "RadmonSubPhysicsListWithLabelMessenger.hh"



RadmonSubPhysicsListWithLabelMessenger *        RadmonSubPhysicsListWithLabelMessenger :: Instance(void)
{
 if (instance==0)
 {
  instance=new RadmonSubPhysicsListWithLabelMessenger;
  
  if (!instance)
   G4Exception("RadmonSubPhysicsListWithLabelMessenger::Instance: RadmonSubPhysicsListWithLabelMessenger singleton not allocated.");
 }    
 
 return instance;
}




  
void                                            RadmonSubPhysicsListWithLabelMessenger :: AddAvailablePhysicsList(const G4String & name)
{
 availablePhysicsLists.insert(name);
}



void                                            RadmonSubPhysicsListWithLabelMessenger :: RemoveAvailablePhysicsList(const G4String & name)
{
 availablePhysicsLists.erase(name);
}





                                                RadmonSubPhysicsListWithLabelMessenger :: RadmonSubPhysicsListWithLabelMessenger()
:
 RadmonMessenger(COMMANDS_PATH, "Interactive labelled physics lists factory commands."),
 RADMON_INITIALIZE_COMMAND(Dump)
{
 RADMON_CREATE_COMMAND_0ARGS(Dump, "Print out the current available physics lists");
}



                                                RadmonSubPhysicsListWithLabelMessenger :: ~RadmonSubPhysicsListWithLabelMessenger()
{
 RADMON_DESTROY_COMMAND(Dump);
}





G4String                                        RadmonSubPhysicsListWithLabelMessenger :: GetCurrentValue(G4UIcommand * /* command */)
{
 G4cout << "RadmonSubPhysicsListWithLabelMessenger::GetCurrentValue(): Not supported" << G4endl;
 
 return G4String();
}



void                                            RadmonSubPhysicsListWithLabelMessenger :: SetNewValue(G4UIcommand * command, G4String newValue)
{
 RADMON_BEGIN_LIST_SET_COMMANDS
  RADMON_SET_COMMAND(Dump)
 RADMON_END_LIST_SET_COMMANDS
}





// Events
void                                            RadmonSubPhysicsListWithLabelMessenger :: OnDump(const G4String & /* value */)
{
 G4String indent("- ");
 
 G4String indent2(indent);
 indent2.prepend("  ");
 
 G4cout << indent << "Available physics list:\n";
 
 AvailablePhysicsLists::const_iterator i(availablePhysicsLists.begin());
 AvailablePhysicsLists::const_iterator end(availablePhysicsLists.end());
 
 if (i==end)
  G4cout << indent2 << "No physics list available.\n";
 else
 {
  while (i!=end)
  {
   G4cout << indent2 << (*i) << '\n';
   i++;
  }
 }
 
 G4cout << G4endl;
}





RadmonSubPhysicsListWithLabelMessenger *        RadmonSubPhysicsListWithLabelMessenger :: instance(0);
