//
// File name:     RadmonSubPhysicsListWithLabelMessenger.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonSubPhysicsListWithLabelMessenger.hh,v 1.2 2005-11-10 08:14:10 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   UI commands for dumping the available constructors
//

#ifndef   RADMONSUBPHYSICSLISTWITHLABELMESSENGER_HH
 #define  RADMONSUBPHYSICSLISTWITHLABELMESSENGER_HH

 // Include files
 #include "RadmonMessenger.hh"
 #include "G4String.hh"
 #include <set>

 class RadmonSubPhysicsListWithLabelMessenger : public RadmonMessenger
 {
  public:
   static RadmonSubPhysicsListWithLabelMessenger * Instance(void);
  
   void                                         AddAvailablePhysicsList(const G4String & name);
   void                                         RemoveAvailablePhysicsList(const G4String & name);
  
   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Hidden constructors and operators
                                                RadmonSubPhysicsListWithLabelMessenger();
                                                RadmonSubPhysicsListWithLabelMessenger(const RadmonSubPhysicsListWithLabelMessenger & copy);
                                               ~RadmonSubPhysicsListWithLabelMessenger();
   RadmonSubPhysicsListWithLabelMessenger &     operator=(const RadmonSubPhysicsListWithLabelMessenger & copy);

  // Private Data Types
   typedef std::set<G4String>                   AvailablePhysicsLists;
   
  // Private variables
   AvailablePhysicsLists                        availablePhysicsLists;
   
   static RadmonSubPhysicsListWithLabelMessenger * instance;
   
  // Commands
   RADMON_DECLARE_COMMAND(Dump);
 };
#endif /* RADMONSUBPHYSICSLISTWITHLABELMESSENGER_HH */
