//
// File name:     RadmonDetectorMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMessenger.hh,v 1.3 2005-09-19 19:42:13 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   UI commands for managing a layout
//

#ifndef   RADMONDETECTORMESSENGER_HH
 #define  RADMONDETECTORMESSENGER_HH

 // Include files
 #include "RadmonMessenger.hh"
 #include "G4String.hh"

 // Forward declarations
 class RadmonVDetectorLayout;
 
 class RadmonDetectorMessenger : public RadmonMessenger
 {
  public:
                                                RadmonDetectorMessenger(RadmonVDetectorLayout * layout);
                                               ~RadmonDetectorMessenger();

   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Hidden constructors and operators
                                                RadmonDetectorMessenger();
                                                RadmonDetectorMessenger(const RadmonDetectorMessenger & copy);
   RadmonDetectorMessenger &                    operator=(const RadmonDetectorMessenger & copy);

  // Private Attributes
   RadmonVDetectorLayout *                      detectorLayout;

  // Commands
   RADMON_DECLARE_COMMAND(EnableEnvironment);
   RADMON_DECLARE_COMMAND(DisableEnvironment);
   RADMON_DECLARE_COMMAND(SetEnvironmentType);
   RADMON_DECLARE_COMMAND(SetEnvironmentAttribute);
   RADMON_DECLARE_COMMAND(ClearEnvironmentAttribute);

   RADMON_DECLARE_COMMAND(CreateMultilayer);
   RADMON_DECLARE_COMMAND(RemoveMultilayer);
   RADMON_DECLARE_COMMAND(SetMultilayerWidth);
   RADMON_DECLARE_COMMAND(SetMultilayerHeight);
   RADMON_DECLARE_COMMAND(AppendLayerToMultilayer);
   RADMON_DECLARE_COMMAND(RemoveLayerFromMultilayer);
   RADMON_DECLARE_COMMAND(RemoveAllLayersFromMultilayer);
   
   RADMON_DECLARE_COMMAND(SetLayerThickness);
   RADMON_DECLARE_COMMAND(SetLayerType);
   RADMON_DECLARE_COMMAND(SetLayerAttribute);
   RADMON_DECLARE_COMMAND(ClearLayerAttribute);

   RADMON_DECLARE_COMMAND(CreatePlacement);
   RADMON_DECLARE_COMMAND(RemovePlacement);
   RADMON_DECLARE_COMMAND(SetPlacementPosition);
   RADMON_DECLARE_COMMAND(SetPlacementRotation);
   RADMON_DECLARE_COMMAND(SetRelativePlacementPosition);
   RADMON_DECLARE_COMMAND(SetRelativePlacementRotation);
   
   RADMON_DECLARE_COMMAND(DumpLayout);
   RADMON_DECLARE_COMMAND(Load);
   RADMON_DECLARE_COMMAND(Save);
 };
#endif /* RADMONDETECTORMESSENGER_HH */
