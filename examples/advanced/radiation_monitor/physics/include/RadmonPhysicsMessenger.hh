//
// File name:     RadmonPhysicsMessenger.hh
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsMessenger.hh,v 1.1 2005-11-07 17:53:48 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   UI commands for managing a layout
//

#ifndef   RADMONPHYSICSMESSENGER_HH
 #define  RADMONPHYSICSMESSENGER_HH

 // Include files
 #include "RadmonMessenger.hh"
 #include "G4String.hh"

 // Forward declarations
 class RadmonVPhysicsLayout;
 
 class RadmonPhysicsMessenger : public RadmonMessenger
 {
  public:
                                                RadmonPhysicsMessenger(RadmonVPhysicsLayout * layout);
                                               ~RadmonPhysicsMessenger();

   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Hidden constructors and operators
                                                RadmonPhysicsMessenger();
                                                RadmonPhysicsMessenger(const RadmonPhysicsMessenger & copy);
   RadmonPhysicsMessenger &                     operator=(const RadmonPhysicsMessenger & copy);

  // Private Attributes
   RadmonVPhysicsLayout *                       physicsLayout;

  // Commands
   RADMON_DECLARE_COMMAND(AddPhysicsList);
   RADMON_DECLARE_COMMAND(RemovePhysicsList);

   RADMON_DECLARE_COMMAND(SetPhysicsListAttribute);
   RADMON_DECLARE_COMMAND(ClearPhysicsListAttribute);

   RADMON_DECLARE_COMMAND(DumpLayout);
   RADMON_DECLARE_COMMAND(Load);
   RADMON_DECLARE_COMMAND(Save);
 };
#endif /* RADMONPHYSICSMESSENGER_HH */
