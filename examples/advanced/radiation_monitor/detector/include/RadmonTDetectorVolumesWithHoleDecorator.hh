//
// File name:     RadmonTDetectorVolumesWithHoleDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorVolumesWithHoleDecorator.hh,v 1.3 2005-09-27 13:54:12 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Decorates with a hole a component with the following method
//                 - Constructor
//                 - Destructor
//                 - RadmonDetectorLayerVolumesList * GenerateVolumesList(void)
//

#ifndef   RADMONTDETECTORVOLUMESWITHHOLEDECORATOR_HH
 #define  RADMONTDETECTORVOLUMESWITHHOLEDECORATOR_HH

 // Include files
 #include "RadmonDetectorLayerVolumeItemSubtraction.hh"

 // Forward declaration
 class RadmonVDetectorLabelledEntityConstructor;
 class RadmonDetectorLayerVolumesList;
 class G4VSolid;

 template <class LayerVolumesComponent>
 class RadmonTDetectorVolumesWithHoleDecorator
 {
  public:
                                                RadmonTDetectorVolumesWithHoleDecorator(const RadmonVDetectorLabelledEntityConstructor * constructor);
                                               ~RadmonTDetectorVolumesWithHoleDecorator();

   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorVolumesWithHoleDecorator();
                                                RadmonTDetectorVolumesWithHoleDecorator(const  RadmonTDetectorVolumesWithHoleDecorator & copy);
   RadmonTDetectorVolumesWithHoleDecorator &    operator=(const  RadmonTDetectorVolumesWithHoleDecorator & copy);

  // Private attributes
   const RadmonVDetectorLabelledEntityConstructor * owner;
   LayerVolumesComponent                        component;
   RadmonDetectorLayerVolumeItemSubtraction     operation;
 };
 
 // Inline implementations
 #include "RadmonTDetectorVolumesWithHoleDecorator.icc"
#endif /* RADMONTDETECTORVOLUMESWITHHOLEDECORATOR_HH */
