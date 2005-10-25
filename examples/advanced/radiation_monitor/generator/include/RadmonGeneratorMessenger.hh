//
// File name:     RadmonGeneratorLayout.hh
// Creation date: Oct 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonGeneratorMessenger.hh,v 1.1 2005-10-25 16:36:41 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Messenger to change sources properties
//

#ifndef   RADMONGENERATORMESSENGER_HH
 #define  RADMONGENERATORMESSENGER_HH
 
 // Include files
 #include "RadmonMessenger.hh"
 
 // Forward declarations
 class RadmonVGeneratorLayout;
 
 class RadmonGeneratorMessenger : public RadmonMessenger
 {
  public:
                                                RadmonGeneratorMessenger(RadmonVGeneratorLayout * layout);
                                               ~RadmonGeneratorMessenger();

   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  // Hidden constructors and operators
  private:
                                                RadmonGeneratorMessenger(void);
                                                RadmonGeneratorMessenger(const RadmonGeneratorMessenger & copy);
   RadmonGeneratorMessenger &                   operator=(const RadmonGeneratorMessenger & copy);

  // Private Attributes
   RadmonVGeneratorLayout *                     generatorLayout;

  // Commands
   RADMON_DECLARE_COMMAND(InsertSource);
   RADMON_DECLARE_COMMAND(SetRelativeSourceIntensity);
   RADMON_DECLARE_COMMAND(RemoveSource);

   RADMON_DECLARE_COMMAND(AppendSourceAlgorithm);
   RADMON_DECLARE_COMMAND(SetSourceAlgorithmType);
   RADMON_DECLARE_COMMAND(RemoveSourceAlgorithm);

   RADMON_DECLARE_COMMAND(SetSourceAlgorithmAttribute);
   RADMON_DECLARE_COMMAND(ClearSourceAlgorithmAttribute);

   RADMON_DECLARE_COMMAND(Load);
   RADMON_DECLARE_COMMAND(Save);
   RADMON_DECLARE_COMMAND(DumpLayout);
 };
#endif /* RADMONGENERATORMESSENGER_HH */
