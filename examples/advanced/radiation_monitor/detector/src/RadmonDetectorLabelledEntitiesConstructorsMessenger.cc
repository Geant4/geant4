//
// File name:     RadmonDetectorLabelledEntitiesConstructorsMessenger.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLabelledEntitiesConstructorsMessenger.cc,v 1.3 2005-11-24 02:35:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Messenger commands path
#define COMMANDS_PATH "/radmon/detectorFactory/"

// Include files
#include "RadmonDetectorLabelledEntitiesConstructorsMessenger.hh"



RadmonDetectorLabelledEntitiesConstructorsMessenger * RadmonDetectorLabelledEntitiesConstructorsMessenger :: Instance(void)
{
 if (instance==0)
 {
  instance=new RadmonDetectorLabelledEntitiesConstructorsMessenger;
  
  if (!instance)
   G4Exception("RadmonDetectorLabelledEntitiesConstructorsMessenger::Instance: RadmonDetectorLabelledEntitiesConstructorsMessenger singleton not allocated.");
 }    
 
 return instance;
}




  
void                                            RadmonDetectorLabelledEntitiesConstructorsMessenger :: AddAvailableConstructor(const G4String & name)
{
 availableConstructors.insert(name);
}



void                                            RadmonDetectorLabelledEntitiesConstructorsMessenger :: RemoveAvailableConstructor(const G4String & name)
{
 availableConstructors.erase(name);
}





                                                RadmonDetectorLabelledEntitiesConstructorsMessenger :: RadmonDetectorLabelledEntitiesConstructorsMessenger()
:
 RadmonMessenger(COMMANDS_PATH, "Interactive labelled detector entities constructors factory commands."),
 RADMON_INITIALIZE_COMMAND(Dump)
{
 RADMON_CREATE_COMMAND_0ARGS(Dump, "Print out the current available detector entity constructors");
}



                                                RadmonDetectorLabelledEntitiesConstructorsMessenger :: ~RadmonDetectorLabelledEntitiesConstructorsMessenger()
{
 RADMON_DESTROY_COMMAND(Dump);
}





G4String                                        RadmonDetectorLabelledEntitiesConstructorsMessenger :: GetCurrentValue(G4UIcommand * /* command */)
{
 G4cout << "RadmonDetectorLabelledEntitiesConstructorsMessenger::GetCurrentValue(): Not supported" << G4endl;
 
 return G4String();
}



void                                            RadmonDetectorLabelledEntitiesConstructorsMessenger :: SetNewValue(G4UIcommand * command, G4String newValue)
{
 RADMON_BEGIN_LIST_SET_COMMANDS
  RADMON_SET_COMMAND(Dump)
 RADMON_END_LIST_SET_COMMANDS
}





// Events
void                                            RadmonDetectorLabelledEntitiesConstructorsMessenger :: OnDump(const G4String & /* value */)
{
 G4String indent("- ");
 
 G4String indent2(indent);
 indent2.prepend("  ");
 
 G4cout << indent << "Available constructors:\n";
 
 AvailableConstructors::const_iterator i(availableConstructors.begin());
 AvailableConstructors::const_iterator end(availableConstructors.end());
 
 if (i==end)
  G4cout << indent2 << "No constructors available.\n";
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





RadmonDetectorLabelledEntitiesConstructorsMessenger * RadmonDetectorLabelledEntitiesConstructorsMessenger :: instance(0);
