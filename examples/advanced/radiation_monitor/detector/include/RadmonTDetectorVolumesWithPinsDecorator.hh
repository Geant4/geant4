//
// File name:     RadmonTDetectorVolumesWithPinsDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorVolumesWithPinsDecorator.hh,v 1.1 2005-09-27 13:55:51 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Decorates with the pin out
//                 - Constructor
//                 - Destructor
//                 - RadmonDetectorLayerVolumesList * GenerateVolumesList(void)
//

#ifndef   RADMONTDETECTORVOLUMESWITHPINSDECORATOR_HH
 #define  RADMONTDETECTORVOLUMESWITHPINSDECORATOR_HH

 // Include files
 #include "RadmonDetectorLayerVolumeItemIntersection.hh"

 // Forward declaration
 class RadmonVDetectorLabelledEntityConstructor;
 class RadmonDetectorLayerVolumesList;

 template <class LayerVolumesComponent>
 class RadmonTDetectorVolumesWithPinsDecorator
 {
  public:
                                                RadmonTDetectorVolumesWithPinsDecorator(const RadmonVDetectorLabelledEntityConstructor * constructor);
                                               ~RadmonTDetectorVolumesWithPinsDecorator();

   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorVolumesWithPinsDecorator();
                                                RadmonTDetectorVolumesWithPinsDecorator(const  RadmonTDetectorVolumesWithPinsDecorator & copy);
   RadmonTDetectorVolumesWithPinsDecorator &    operator=(const  RadmonTDetectorVolumesWithPinsDecorator & copy);

  // Private attributes
   const RadmonVDetectorLabelledEntityConstructor * owner;
   LayerVolumesComponent                        component;
   G4VisAttributes *                            visAttributes;
   RadmonDetectorLayerVolumeItemIntersection    intersection;
 };
 
 // Inline implementations
 #include "RadmonTDetectorVolumesWithPinsDecorator.icc"
#endif /* RADMONTDETECTORVOLUMESWITHPINSDECORATOR_HH */
