//
// File name:     RadmonSubPhysicsListWithLabelMessenger.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSubPhysicsListWithLabelMessenger.cc,v 1.1 2005-11-07 17:52:36 capra Exp $
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
 RadmonMessenger(COMMANDS_PATH, "Interactive labelled detector entities constructors factory commands."),
 RADMON_INITIALIZE_COMMAND(Dump)
{
 RADMON_CREATE_COMMAND_0ARGS(Dump, "Print out the current available detector enetity constructors");
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
