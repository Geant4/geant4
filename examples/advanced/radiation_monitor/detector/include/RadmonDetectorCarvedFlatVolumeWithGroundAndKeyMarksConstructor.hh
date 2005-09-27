//
// File name:     RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor.hh,v 1.1 2005-09-27 13:55:51 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a carved box with a ground pad and marks
//

#ifndef   RADMONDETECTORCARVEDFLATVOLUMEWITHGROUNDANDKEYMARKSCONSTRUCTOR_HH
 #define  RADMONDETECTORCARVEDFLATVOLUMEWITHGROUNDANDKEYMARKSCONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorCarvedVolumesDecorator.hh"
 #include "RadmonTDetectorVolumesWithGroundDecorator.hh"
 #include "RadmonTDetectorVolumesWithKeyMarksDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorCarvedVolumesDecorator<RadmonTDetectorVolumesWithKeyMarksDecorator<RadmonTDetectorVolumesWithGroundDecorator<RadmonDetectorFlatVolumeComponent> > > >
 {
  public:
   inline                                       RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor();
   inline virtual                              ~RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor(const RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor & copy);
   RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor & operator=(const RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor.icc"
#endif /* RADMONDETECTORCARVEDFLATVOLUMEWITHGROUNDANDKEYMARKSCONSTRUCTOR_HH */
