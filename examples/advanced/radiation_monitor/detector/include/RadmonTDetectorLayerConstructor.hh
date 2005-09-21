//
// File name:     RadmonTDetectorLayerConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorLayerConstructor.hh,v 1.2 2005-09-21 14:48:18 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Generates a layer using a component with the following methods
//                 - Constructor with RadmonVDetectorLabelledEntityConstructor or its base classes as argument
//                 - Destructor
//                 - RadmonDetectorLayerVolumesList * GenerateVolumesList(void)
//

#ifndef   RADMONTDETECTORLAYERCONSTRUCTOR_HH
 #define  RADMONTDETECTORLAYERCONSTRUCTOR_HH

 // Include files
 #include "RadmonVDetectorLabelledEntityConstructor.hh"

 #include "globals.hh" 

 // Forward declaration
 class RadmonDetectorLayerVolumesList;

 template <class LayerVolumesComponent>
 class RadmonTDetectorLayerConstructor : public RadmonVDetectorLabelledEntityConstructor
 {
  public:
                                                RadmonTDetectorLayerConstructor(const G4String & label);
                                               ~RadmonTDetectorLayerConstructor();

   virtual G4LogicalVolume *                    ConstructLogicalVolume(void);
   virtual RadmonVDetectorLabelledEntityConstructor * New(void) const;

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorLayerConstructor();
                                                RadmonTDetectorLayerConstructor(const RadmonTDetectorLayerConstructor & copy);
   RadmonTDetectorLayerConstructor &            operator=(const RadmonTDetectorLayerConstructor & copy);

  // Private attributes
   RadmonDetectorLayerVolumesList *             volumesList;
   LayerVolumesComponent                        component;
 };
 
 // Inline implementations
 #include "RadmonTDetectorLayerConstructor.icc"
#endif /* RADMONTDETECTORLAYERCONSTRUCTOR_HH */
