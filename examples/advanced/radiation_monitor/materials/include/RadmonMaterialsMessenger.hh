//
// File name:     RadmonMaterialsMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonMaterialsMessenger.hh,v 1.2 2005-09-19 19:40:26 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   UI commands for managing a layout
//

#ifndef   RADMONMATERIALSMESSENGER_HH
 #define  RADMONMATERIALSMESSENGER_HH

 // Include files
 #include "RadmonMessenger.hh"
 #include "G4String.hh"

 // Forward declarations
 class RadmonMaterialsManager;
 
 class RadmonMaterialsMessenger : public RadmonMessenger
 {
  public:
                                                RadmonMaterialsMessenger(RadmonMaterialsManager * manager);
                                               ~RadmonMaterialsMessenger();

   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Hidden constructors and operators
                                                RadmonMaterialsMessenger();
                                                RadmonMaterialsMessenger(const RadmonMaterialsMessenger & copy);
   RadmonMaterialsMessenger &                   operator=(const RadmonMaterialsMessenger & copy);

  // Private Attributes
   RadmonMaterialsManager *                     materialsManager;

  // Commands
   RADMON_DECLARE_COMMAND(CreateElement);
   RADMON_DECLARE_COMMAND(CreateMaterial);
   RADMON_DECLARE_COMMAND(AddComponentByAtoms);
   RADMON_DECLARE_COMMAND(AddComponentByFraction);
   RADMON_DECLARE_COMMAND(SetMaterialColor);
   RADMON_DECLARE_COMMAND(SetMaterialVisibility);
   RADMON_DECLARE_COMMAND(SetMaterialStyle);
   RADMON_DECLARE_COMMAND(Dump);
   RADMON_DECLARE_COMMAND(Insert);
   RADMON_DECLARE_COMMAND(Save);
 };
#endif /* RADMONMATERIALSMESSENGER_HH */
