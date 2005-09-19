//
// File name:     RadmonDetectorLabelledEntitiesConstructorsMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLabelledEntitiesConstructorsMessenger.hh,v 1.1 2005-09-19 19:38:50 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   UI commands for dumping the available constructors
//

#ifndef   RADMONDETECTORLABELLEDENTITIESCONSTRUCTORSMESSENGER_HH
 #define  RADMONDETECTORLABELLEDENTITIESCONSTRUCTORSMESSENGER_HH

 // Include files
 #include "RadmonMessenger.hh"
 #include "G4String.hh"
 #include <set>

 class RadmonDetectorLabelledEntitiesConstructorsMessenger : public RadmonMessenger
 {
  public:
   static RadmonDetectorLabelledEntitiesConstructorsMessenger * Instance(void);
  
   void                                         AddAvailableConstructor(const G4String & name);
   void                                         RemoveAvailableConstructor(const G4String & name);
  
   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Hidden constructors and operators
                                                RadmonDetectorLabelledEntitiesConstructorsMessenger();
                                                RadmonDetectorLabelledEntitiesConstructorsMessenger(const RadmonDetectorLabelledEntitiesConstructorsMessenger & copy);
                                               ~RadmonDetectorLabelledEntitiesConstructorsMessenger();
   RadmonDetectorLabelledEntitiesConstructorsMessenger & operator=(const RadmonDetectorLabelledEntitiesConstructorsMessenger & copy);

  // Private Data Types
   typedef std::set<G4String>                   AvailableConstructors;
   
  // Private variables
   AvailableConstructors                        availableConstructors;
   
   static RadmonDetectorLabelledEntitiesConstructorsMessenger * instance;
   
  // Commands
   RADMON_DECLARE_COMMAND(Dump);
 };
#endif /* RADMONDETECTORLABELLEDENTITIESCONSTRUCTORSMESSENGER_HH */
