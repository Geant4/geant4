//
// File name:     RadmonDetectorFlatVolumeWithGroundAndKeyMarksConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeWithGroundAndKeyMarksConstructor.hh,v 1.1 2005-09-27 13:55:51 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with a ground pad and marks
//

#ifndef   RADMONDETECTORFLATVOLUMEWITHGROUNDANDKEYMARKSCONSTRUCTOR_HH
 #define  RADMONDETECTORFLATVOLUMEWITHGROUNDANDKEYMARKSCONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorVolumesWithGroundDecorator.hh"
 #include "RadmonTDetectorVolumesWithKeyMarksDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorFlatVolumeWithGroundAndKeyMarksConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorVolumesWithKeyMarksDecorator<RadmonTDetectorVolumesWithGroundDecorator<RadmonDetectorFlatVolumeComponent> > >
 {
  public:
   inline                                       RadmonDetectorFlatVolumeWithGroundAndKeyMarksConstructor();
   inline virtual                              ~RadmonDetectorFlatVolumeWithGroundAndKeyMarksConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorFlatVolumeWithGroundAndKeyMarksConstructor(const RadmonDetectorFlatVolumeWithGroundAndKeyMarksConstructor & copy);
   RadmonDetectorFlatVolumeWithGroundAndKeyMarksConstructor & operator=(const RadmonDetectorFlatVolumeWithGroundAndKeyMarksConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorFlatVolumeWithGroundAndKeyMarksConstructor.icc"
#endif /* RADMONDETECTORFLATVOLUMEWITHGROUNDANDKEYMARKSCONSTRUCTOR_HH */
