//
// File name:     RadmonDetectorFlatVolumeWithPadsConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeWithPadsConstructor.hh,v 1.1 2005-09-30 08:33:30 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with pins
//

#ifndef   RADMONDETECTORFLATVOLUMEWITHPADSCONSTRUCTOR_HH
 #define  RADMONDETECTORFLATVOLUMEWITHPADSCONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorVolumesWithPadsDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorFlatVolumeWithPadsConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorVolumesWithPadsDecorator<RadmonDetectorFlatVolumeComponent> >
 {
  public:
   inline                                       RadmonDetectorFlatVolumeWithPadsConstructor();
   inline virtual                              ~RadmonDetectorFlatVolumeWithPadsConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorFlatVolumeWithPadsConstructor(const RadmonDetectorFlatVolumeWithPadsConstructor & copy);
   RadmonDetectorFlatVolumeWithPadsConstructor & operator=(const RadmonDetectorFlatVolumeWithPadsConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorFlatVolumeWithPadsConstructor.icc"
#endif /* RADMONDETECTORFLATVOLUMEWITHPADSCONSTRUCTOR_HH */
