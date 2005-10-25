//
// File name:     RadmonGeneratorsWithLabelMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorsWithLabelMessenger.hh,v 1.1 2005-10-25 16:36:41 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   UI commands for dumping the available constructors
//

#ifndef   RADMONGENERATORSWITHLABELMESSENGER_HH
 #define  RADMONGENERATORSWITHLABELMESSENGER_HH

 // Include files
 #include "RadmonMessenger.hh"
 #include "G4String.hh"
 #include <set>

 class RadmonGeneratorsWithLabelMessenger : public RadmonMessenger
 {
  public:
   static RadmonGeneratorsWithLabelMessenger *  Instance(void);
  
   void                                         AddAvailableGenerator(const G4String & name);
   void                                         RemoveAvailableGenerator(const G4String & name);
  
   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Hidden constructors and operators
                                                RadmonGeneratorsWithLabelMessenger();
                                                RadmonGeneratorsWithLabelMessenger(const RadmonGeneratorsWithLabelMessenger & copy);
                                               ~RadmonGeneratorsWithLabelMessenger();
   RadmonGeneratorsWithLabelMessenger &         operator=(const RadmonGeneratorsWithLabelMessenger & copy);

  // Private Data Types
   typedef std::set<G4String>                   AvailableGenerators;
   
  // Private variables
   AvailableGenerators                          availableGenerators;
   
   static RadmonGeneratorsWithLabelMessenger * instance;
   
  // Commands
   RADMON_DECLARE_COMMAND(Dump);
 };
#endif /* RADMONGENERATORSWITHLABELMESSENGER_HH */
