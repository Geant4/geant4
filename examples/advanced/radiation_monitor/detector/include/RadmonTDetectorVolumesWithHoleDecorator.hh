//
// File name:     RadmonTDetectorVolumesWithHoleDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorVolumesWithHoleDecorator.hh,v 1.2 2005-09-21 14:52:57 capra Exp $
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
 #include <stack>

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

  // Private data types
   typedef std::stack<G4VSolid *>               OwnedSolids;
  
  // Private attributes
   const RadmonVDetectorLabelledEntityConstructor * owner;
   LayerVolumesComponent                        component;
   OwnedSolids                                  ownedSolids;
 };
 
 // Inline implementations
 #include "RadmonTDetectorVolumesWithHoleDecorator.icc"
#endif /* RADMONTDETECTORVOLUMESWITHHOLEDECORATOR_HH */
