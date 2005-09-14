//
// File name:     RadmonDetectorMessenger.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorMessenger.hh,v 1.2 2005-09-14 12:28:31 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   UI commands for managing a layout
//

#ifndef   RADMONDETECTORMESSENGER_HH
 #define  RADMONDETECTORMESSENGER_HH

 // Include files
 #include "G4UImessenger.hh"
 #include "G4String.hh"

 // Forward declarations
 class RadmonVDetectorLayout;
 class G4UIdirectory;
 class G4UIcommand;
 
 #define RADMONDETECTORMESSENGER_DECLARE_COMMAND(command) private:                                      \
                                                           void On ## command (const G4String & value); \
                                                           G4UIcommand * cmd ## command
 
 class RadmonDetectorMessenger : public G4UImessenger
 {
  public:
                                                RadmonDetectorMessenger(RadmonVDetectorLayout * layout);
                                               ~RadmonDetectorMessenger();

   virtual G4String                             GetCurrentValue(G4UIcommand * command);
   virtual void                                 SetNewValue(G4UIcommand * command, G4String newValue);

  private:
  // Private methods
   G4bool                                       ProcessArguments(const G4String & rawArguments, G4int nArgs, G4String * arguments) const;
   G4double                                     GetUnit(const G4String & unitStr, const char *cathegory) const;
  
  // Hidden constructors and operators
                                                RadmonDetectorMessenger();
                                                RadmonDetectorMessenger(const RadmonDetectorMessenger & copy);
   RadmonDetectorMessenger &                    operator=(const RadmonDetectorMessenger & copy);

  // Private Attributes
   RadmonVDetectorLayout *                      detectorLayout;

   G4UIdirectory *                              directory;

  // Commands
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(EnableEnvironment);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(DisableEnvironment);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(SetEnvironmentType);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(SetEnvironmentAttribute);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(ClearEnvironmentAttribute);

   RADMONDETECTORMESSENGER_DECLARE_COMMAND(CreateMultilayer);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(RemoveMultilayer);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(SetMultilayerWidth);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(SetMultilayerHeight);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(AppendLayerToMultilayer);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(RemoveLayerFromMultilayer);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(RemoveAllLayersFromMultilayer);
   
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(SetLayerThickness);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(SetLayerType);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(SetLayerAttribute);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(ClearLayerAttribute);

   RADMONDETECTORMESSENGER_DECLARE_COMMAND(CreatePlacement);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(RemovePlacement);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(SetPlacementPosition);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(SetPlacementRotation);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(SetRelativePlacementPosition);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(SetRelativePlacementRotation);
   
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(DumpLayout);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(Load);
   RADMONDETECTORMESSENGER_DECLARE_COMMAND(Save);
 };
#endif /* RADMONDETECTORMESSENGER_HH */
