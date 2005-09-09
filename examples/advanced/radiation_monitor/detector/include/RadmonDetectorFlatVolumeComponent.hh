//
// File name:     RadmonDetectorFlatVolumeComponent.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeComponent.hh,v 1.1 2005-09-09 08:26:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Component to create a uniform layer
//

#ifndef   RADMONDETECTORFLATVOLUMECOMPONENT_HH
 #define  RADMONDETECTORFLATVOLUMECOMPONENT_HH
 
 // Include files
 #include "globals.hh"
 
 // Forward declarations
 class RadmonTDetectorLayerConstructor;
 class RadmonDetectorLayerVolumesList;
 class G4box;
 
 class RadmonDetectorFlatVolumeComponent
 {
  public:
                                                RadmonDetectorFlatVolumeComponent(const RadmonTDetectorLayerConstructor * owner);
                                               ~RadmonDetectorFlatVolumeComponent();
   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonDetectorFlatVolumeComponent();
                                                RadmonDetectorFlatVolumeComponent(const RadmonDetectorFlatVolumeComponent & copy);
   RadmonDetectorFlatVolumeComponent &          operator=(const RadmonDetectorFlatVolumeComponent & copy);

  // Private attributes
   G4double                                     width;
   G4double                                     height;
   G4double                                     thickness;
   G4Box *                                      box;
 };
#endif /* RADMONDETECTORFLATVOLUMECOMPONENT_HH */
