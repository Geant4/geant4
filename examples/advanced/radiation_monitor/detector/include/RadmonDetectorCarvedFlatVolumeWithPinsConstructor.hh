//
// File name:     RadmonDetectorCarvedFlatVolumeWithPinsConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorCarvedFlatVolumeWithPinsConstructor.hh,v 1.1 2005-09-27 13:55:51 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with carved borders and pins
//

#ifndef   RADMONDETECTORCARVEDFLATVOLUMEWITHPINSCONSTRUCTOR_HH
 #define  RADMONDETECTORCARVEDFLATVOLUMEWITHPINSCONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorCarvedVolumesDecorator.hh"
 #include "RadmonTDetectorVolumesWithPinsDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorCarvedFlatVolumeWithPinsConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorCarvedVolumesDecorator<RadmonTDetectorVolumesWithPinsDecorator<RadmonDetectorFlatVolumeComponent> > >
 {
  public:
   inline                                       RadmonDetectorCarvedFlatVolumeWithPinsConstructor();
   inline virtual                              ~RadmonDetectorCarvedFlatVolumeWithPinsConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorCarvedFlatVolumeWithPinsConstructor(const RadmonDetectorCarvedFlatVolumeWithPinsConstructor & copy);
   RadmonDetectorCarvedFlatVolumeWithPinsConstructor & operator=(const RadmonDetectorCarvedFlatVolumeWithPinsConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorCarvedFlatVolumeWithPinsConstructor.icc"
#endif /* RADMONDETECTORCARVEDFLATVOLUMEWITHPINSCONSTRUCTOR_HH */
