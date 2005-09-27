//
// File name:     RadmonDetectorFlatVolumeWithGroundConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeWithGroundConstructor.hh,v 1.1 2005-09-27 13:55:51 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with carved borders
//

#ifndef   RADMONDETECTORFLATVOLUMEWITHGROUNDCONSTRUCTOR_HH
 #define  RADMONDETECTORFLATVOLUMEWITHGROUNDCONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorVolumesWithGroundDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorFlatVolumeWithGroundConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorVolumesWithGroundDecorator<RadmonDetectorFlatVolumeComponent> >
 {
  public:
   inline                                       RadmonDetectorFlatVolumeWithGroundConstructor();
   inline virtual                              ~RadmonDetectorFlatVolumeWithGroundConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorFlatVolumeWithGroundConstructor(const RadmonDetectorFlatVolumeWithGroundConstructor & copy);
   RadmonDetectorFlatVolumeWithGroundConstructor & operator=(const RadmonDetectorFlatVolumeWithGroundConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorFlatVolumeWithGroundConstructor.icc"
#endif /* RADMONDETECTORFLATVOLUMEWITHGROUNDCONSTRUCTOR_HH */
