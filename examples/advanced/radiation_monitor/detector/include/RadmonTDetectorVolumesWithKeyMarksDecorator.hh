//
// File name:     RadmonTDetectorVolumesWithKeyMarksDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorVolumesWithKeyMarksDecorator.hh,v 1.1 2005-09-27 13:55:51 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Decorates with key marks
//                 - Constructor
//                 - Destructor
//                 - RadmonDetectorLayerVolumesList * GenerateVolumesList(void)
//

#ifndef   RADMONTDETECTORVOLUMESWITHKEYMARKSDECORATOR_HH
 #define  RADMONTDETECTORVOLUMESWITHKEYMARKSDECORATOR_HH

 // Include files
 #include "RadmonDetectorLayerVolumeItemIntersection.hh"

 // Forward declaration
 class RadmonVDetectorLabelledEntityConstructor;
 class RadmonDetectorLayerVolumesList;
 class G4VSolid;

 template <class LayerVolumesComponent>
 class RadmonTDetectorVolumesWithKeyMarksDecorator
 {
  public:
                                                RadmonTDetectorVolumesWithKeyMarksDecorator(const RadmonVDetectorLabelledEntityConstructor * constructor);
                                               ~RadmonTDetectorVolumesWithKeyMarksDecorator();

   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorVolumesWithKeyMarksDecorator();
                                                RadmonTDetectorVolumesWithKeyMarksDecorator(const  RadmonTDetectorVolumesWithKeyMarksDecorator & copy);
   RadmonTDetectorVolumesWithKeyMarksDecorator & operator=(const  RadmonTDetectorVolumesWithKeyMarksDecorator & copy);

  // Private attributes
   const RadmonVDetectorLabelledEntityConstructor * owner;
   LayerVolumesComponent                        component;
   G4VisAttributes *                            visAttributes;
   G4VSolid *                                   box;
   G4VSolid *                                   tub;
   G4RotationMatrix                             identity;
   RadmonDetectorLayerVolumeItemIntersection    intersection;
 };
 
 // Inline implementations
 #include "RadmonTDetectorVolumesWithKeyMarksDecorator.icc"
#endif /* RADMONTDETECTORVOLUMESWITHKEYMARKSDECORATOR_HH */
