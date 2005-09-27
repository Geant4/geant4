//
// File name:     RadmonDetectorFlatVolumeWithPinsConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeWithPinsConstructor.hh,v 1.1 2005-09-27 13:55:51 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with pins
//

#ifndef   RADMONDETECTORFLATVOLUMEWITHPINSCONSTRUCTOR_HH
 #define  RADMONDETECTORFLATVOLUMEWITHPINSCONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorVolumesWithPinsDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorFlatVolumeWithPinsConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorVolumesWithPinsDecorator<RadmonDetectorFlatVolumeComponent> >
 {
  public:
   inline                                       RadmonDetectorFlatVolumeWithPinsConstructor();
   inline virtual                              ~RadmonDetectorFlatVolumeWithPinsConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorFlatVolumeWithPinsConstructor(const RadmonDetectorFlatVolumeWithPinsConstructor & copy);
   RadmonDetectorFlatVolumeWithPinsConstructor & operator=(const RadmonDetectorFlatVolumeWithPinsConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorFlatVolumeWithPinsConstructor.icc"
#endif /* RADMONDETECTORFLATVOLUMEWITHPINSCONSTRUCTOR_HH */
