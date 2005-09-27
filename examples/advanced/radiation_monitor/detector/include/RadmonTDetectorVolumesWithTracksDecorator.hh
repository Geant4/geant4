//
// File name:     RadmonTDetectorVolumesWithTracksDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorVolumesWithTracksDecorator.hh,v 1.1 2005-09-27 13:55:51 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Decorates with tracks according to the pin-out
//                 - Constructor
//                 - Destructor
//                 - RadmonDetectorLayerVolumesList * GenerateVolumesList(void)
//

#ifndef   RADMONTDETECTORVOLUMESWITHTRACKSDECORATOR_HH
 #define  RADMONTDETECTORVOLUMESWITHTRACKSDECORATOR_HH

 // Include files
 #include "RadmonDetectorLayerVolumeItemIntersection.hh"

 // Forward declaration
 class RadmonVDetectorLabelledEntityConstructor;
 class RadmonDetectorLayerVolumesList;
 class G4VSolid;

 template <class LayerVolumesComponent>
 class RadmonTDetectorVolumesWithTracksDecorator
 {
  public:
                                                RadmonTDetectorVolumesWithTracksDecorator(const RadmonVDetectorLabelledEntityConstructor * constructor);
                                               ~RadmonTDetectorVolumesWithTracksDecorator();

   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorVolumesWithTracksDecorator();
                                                RadmonTDetectorVolumesWithTracksDecorator(const  RadmonTDetectorVolumesWithTracksDecorator & copy);
   RadmonTDetectorVolumesWithTracksDecorator &  operator=(const  RadmonTDetectorVolumesWithTracksDecorator & copy);

  // Private attributes
   const RadmonVDetectorLabelledEntityConstructor * owner;
   LayerVolumesComponent                        component;
   G4VisAttributes *                            visAttributes;
   RadmonDetectorLayerVolumeItemIntersection    intersection;
 };
 
 // Inline implementations
 #include "RadmonTDetectorVolumesWithTracksDecorator.icc"
#endif /* RADMONTDETECTORVOLUMESWITHTRACKSDECORATOR_HH */
