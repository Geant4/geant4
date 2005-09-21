//
// File name:     RadmonDetectorCarvedFlatVolumeWithHoleConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorCarvedFlatVolumeWithHoleConstructor.hh,v 1.1 2005-09-21 14:50:42 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with carved borders
//

#ifndef   RADMONDETECTORCARVEDFLATVOLUMEWITHHOLECONSTRUCTOR_HH
 #define  RADMONDETECTORCARVEDFLATVOLUMEWITHHOLECONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorCarvedVolumesDecorator.hh"
 #include "RadmonTDetectorVolumesWithHoleDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorCarvedFlatVolumeWithHoleConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorCarvedVolumesDecorator<RadmonTDetectorVolumesWithHoleDecorator<RadmonDetectorFlatVolumeComponent> > >
 {
  public:
   inline                                       RadmonDetectorCarvedFlatVolumeWithHoleConstructor();
   inline virtual                              ~RadmonDetectorCarvedFlatVolumeWithHoleConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorCarvedFlatVolumeWithHoleConstructor(const RadmonDetectorCarvedFlatVolumeWithHoleConstructor & copy);
   RadmonDetectorCarvedFlatVolumeWithHoleConstructor & operator=(const RadmonDetectorCarvedFlatVolumeWithHoleConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorCarvedFlatVolumeWithHoleConstructor.icc"
#endif /* RADMONDETECTORCARVEDFLATVOLUMEWITHHOLECONSTRUCTOR_HH */
