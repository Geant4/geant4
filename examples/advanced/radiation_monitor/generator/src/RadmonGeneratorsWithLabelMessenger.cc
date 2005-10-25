//
// File name:     RadmonGeneratorsWithLabelMessenger.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorsWithLabelMessenger.cc,v 1.1 2005-10-25 16:36:43 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Messenger commands path
#define COMMANDS_PATH "/radmon/generatorsFactory/"

// Include files
#include "RadmonGeneratorsWithLabelMessenger.hh"



RadmonGeneratorsWithLabelMessenger *            RadmonGeneratorsWithLabelMessenger :: Instance(void)
{
 if (instance==0)
 {
  instance=new RadmonGeneratorsWithLabelMessenger;
  
  if (!instance)
   G4Exception("RadmonGeneratorsWithLabelMessenger::Instance: RadmonGeneratorsWithLabelMessenger singleton not allocated.");
 }    
 
 return instance;
}




  
void                                            RadmonGeneratorsWithLabelMessenger :: AddAvailableGenerator(const G4String & name)
{
 availableGenerators.insert(name);
}



void                                            RadmonGeneratorsWithLabelMessenger :: RemoveAvailableGenerator(const G4String & name)
{
 availableGenerators.erase(name);
}





                                                RadmonGeneratorsWithLabelMessenger :: RadmonGeneratorsWithLabelMessenger()
:
 RadmonMessenger(COMMANDS_PATH, "Interactive primary generators factory commands."),
 RADMON_INITIALIZE_COMMAND(Dump)
{
 RADMON_CREATE_COMMAND_0ARGS(Dump, "Print out the current available primary generators");
}



                                                RadmonGeneratorsWithLabelMessenger :: ~RadmonGeneratorsWithLabelMessenger()
{
 RADMON_DESTROY_COMMAND(Dump);
}





G4String                                        RadmonGeneratorsWithLabelMessenger :: GetCurrentValue(G4UIcommand * /* command */)
{
 G4cout << "RadmonGeneratorsWithLabelMessenger::GetCurrentValue(): Not supported" << G4endl;
 
 return G4String();
}



void                                            RadmonGeneratorsWithLabelMessenger :: SetNewValue(G4UIcommand * command, G4String newValue)
{
 RADMON_BEGIN_LIST_SET_COMMANDS
  RADMON_SET_COMMAND(Dump)
 RADMON_END_LIST_SET_COMMANDS
}





// Events
void                                            RadmonGeneratorsWithLabelMessenger :: OnDump(const G4String & /* value */)
{
 G4String indent("- ");
 
 G4String indent2(indent);
 indent2.prepend("  ");
 
 G4cout << indent << "Available generators:\n";
 
 AvailableGenerators::const_iterator i(availableGenerators.begin());
 AvailableGenerators::const_iterator end(availableGenerators.end());
 
 if (i==end)
  G4cout << indent2 << "No generators available.\n";
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





RadmonGeneratorsWithLabelMessenger *            RadmonGeneratorsWithLabelMessenger :: instance(0);
