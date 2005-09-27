//
// File name:     RadmonTDetectorCarvedVolumesDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorCarvedVolumesDecorator.hh,v 1.3 2005-09-27 13:54:12 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Carves the borders of a component with the following method
//                 - Constructor
//                 - Destructor
//                 - RadmonDetectorLayerVolumesList * GenerateVolumesList(void)
//

#ifndef   RADMONTDETECTORCARVEDVOLUMESDECORATOR_HH
 #define  RADMONTDETECTORCARVEDVOLUMESDECORATOR_HH

 // Include files
 #include <stack>
 #include "RadmonDetectorLayerVolumeItemSubtraction.hh"
 #include "RadmonDetectorLayerVolumeItemIntersection.hh"
 #include "G4RotationMatrix.hh"

 // Forward declaration
 class RadmonVDetectorLabelledEntityConstructor;
 class RadmonDetectorLayerVolumesList;
 class G4VSolid;
 class G4VAttributes;

 template <class LayerVolumesComponent>
 class RadmonTDetectorCarvedVolumesDecorator
 {
  public:
                                                RadmonTDetectorCarvedVolumesDecorator(const RadmonVDetectorLabelledEntityConstructor * constructor);
                                               ~RadmonTDetectorCarvedVolumesDecorator();

   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorCarvedVolumesDecorator();
                                                RadmonTDetectorCarvedVolumesDecorator(const RadmonTDetectorCarvedVolumesDecorator & copy);
   RadmonTDetectorCarvedVolumesDecorator &      operator=(const RadmonTDetectorCarvedVolumesDecorator & copy);

  // Private data types
   typedef std::stack<G4VSolid *>               OwnedSolids;
  
  // Private attributes
   const RadmonVDetectorLabelledEntityConstructor * owner;
   LayerVolumesComponent                        component;
   OwnedSolids                                  ownedSolids;
   G4VisAttributes *                            visAttributes;
   RadmonDetectorLayerVolumeItemSubtraction     subtraction;
   RadmonDetectorLayerVolumeItemIntersection    intersection;
   G4RotationMatrix                             identity;
  };
 
 // Inline implementations
 #include "RadmonTDetectorCarvedVolumesDecorator.icc"
#endif /* RADMONTDETECTORCARVEDVOLUMESDECORATOR_HH */
