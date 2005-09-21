//
// File name:     RadmonDetectorFlatVolumeComponent.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeComponent.hh,v 1.2 2005-09-21 14:52:02 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Component to create a box
//

#ifndef   RADMONDETECTORFLATVOLUMECOMPONENT_HH
 #define  RADMONDETECTORFLATVOLUMECOMPONENT_HH
 
 // Include files
 #include "globals.hh"
 
 // Forward declarations
 class RadmonVDetectorLabelledEntityConstructor;
 class RadmonDetectorLayerVolumesList;
 class G4Box;
 class G4VisAttributes;
 
 class RadmonDetectorFlatVolumeComponent
 {
  public:
   inline                                       RadmonDetectorFlatVolumeComponent(const RadmonVDetectorLabelledEntityConstructor * constructor);
                                               ~RadmonDetectorFlatVolumeComponent();
   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonDetectorFlatVolumeComponent();
                                                RadmonDetectorFlatVolumeComponent(const RadmonDetectorFlatVolumeComponent & copy);
   RadmonDetectorFlatVolumeComponent &          operator=(const RadmonDetectorFlatVolumeComponent & copy);

  // Private attributes
   const RadmonVDetectorLabelledEntityConstructor * owner;
   G4Box *                                      box;
   G4VisAttributes *                            visAttributes;
 };
 
 // Inline implementations
 #include "RadmonDetectorFlatVolumeComponent.icc"
#endif /* RADMONDETECTORFLATVOLUMECOMPONENT_HH */
