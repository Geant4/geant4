//
// File name:     RadmonTDetectorVolumesWithPadsDecorator.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonTDetectorVolumesWithPadsDecorator.hh,v 1.1 2005-09-30 08:33:30 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Decorates with a ground plane
//                 - Constructor
//                 - Destructor
//                 - RadmonDetectorLayerVolumesList * GenerateVolumesList(void)
//

#ifndef   RADMONTDETECTORVOLUMESWITHPADSDECORATOR_HH
 #define  RADMONTDETECTORVOLUMESWITHPADSDECORATOR_HH

 // Include files
 #include "RadmonDetectorLayerVolumeItemIntersection.hh"
 #include "RadmonDetectorPadsData.hh"
 #include <stack>

 // Forward declaration
 class RadmonVDetectorLabelledEntityConstructor;
 class RadmonDetectorLayerVolumesList;
 class G4VSolid;

 template <class LayerVolumesComponent>
 class RadmonTDetectorVolumesWithPadsDecorator
 {
  public:
                                                RadmonTDetectorVolumesWithPadsDecorator(const RadmonVDetectorLabelledEntityConstructor * constructor);
                                               ~RadmonTDetectorVolumesWithPadsDecorator();

   RadmonDetectorLayerVolumesList *             GenerateVolumesList(void);

  private:
  // Hidden constructors and operators
                                                RadmonTDetectorVolumesWithPadsDecorator();
                                                RadmonTDetectorVolumesWithPadsDecorator(const  RadmonTDetectorVolumesWithPadsDecorator & copy);
   RadmonTDetectorVolumesWithPadsDecorator &    operator=(const  RadmonTDetectorVolumesWithPadsDecorator & copy);

  // Private datatypes
   typedef std::stack<G4VisAttributes *>        AttributesStack;
   typedef std::list<RadmonDetectorPadsData>    PadsDataList;

  // Private attributes
   const RadmonVDetectorLabelledEntityConstructor * owner;
   LayerVolumesComponent                        component;
   RadmonDetectorLayerVolumeItemIntersection    intersection;
   AttributesStack                              attributesStack;
   PadsDataList                                 padsDataList;
 };
 
 // Inline implementations
 #include "RadmonTDetectorVolumesWithPadsDecorator.icc"
#endif /* RADMONTDETECTORVOLUMESWITHPADSDECORATOR_HH */
