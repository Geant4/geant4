//
// File name:     RadmonDetectorFlatVolumeWithHoleConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeWithHoleConstructor.hh,v 1.1 2005-09-21 14:50:42 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with a hole
//

#ifndef   RADMONDETECTORFLATVOLUMEWITHHOLECONSTRUCTOR_HH
 #define  RADMONDETECTORFLATVOLUMEWITHHOLECONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorVolumesWithHoleDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorFlatVolumeWithHoleConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorVolumesWithHoleDecorator<RadmonDetectorFlatVolumeComponent> >
 {
  public:
   inline                                       RadmonDetectorFlatVolumeWithHoleConstructor();
   inline virtual                              ~RadmonDetectorFlatVolumeWithHoleConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorFlatVolumeWithHoleConstructor(const RadmonDetectorFlatVolumeWithHoleConstructor & copy);
   RadmonDetectorFlatVolumeWithHoleConstructor & operator=(const RadmonDetectorFlatVolumeWithHoleConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorFlatVolumeWithHoleConstructor.icc"
#endif /* RADMONDETECTORFLATVOLUMEWITHHOLECONSTRUCTOR_HH */
