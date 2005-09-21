//
// File name:     RadmonDetectorCarvedFlatVolumeConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorCarvedFlatVolumeConstructor.hh,v 1.1 2005-09-21 14:50:42 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with carved borders
//

#ifndef   RADMONDETECTORCARVEDFLATVOLUMECONSTRUCTOR_HH
 #define  RADMONDETECTORCARVEDFLATVOLUMECONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorCarvedVolumesDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorCarvedFlatVolumeConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorCarvedVolumesDecorator<RadmonDetectorFlatVolumeComponent> >
 {
  public:
   inline                                       RadmonDetectorCarvedFlatVolumeConstructor();
   inline virtual                              ~RadmonDetectorCarvedFlatVolumeConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorCarvedFlatVolumeConstructor(const RadmonDetectorCarvedFlatVolumeConstructor & copy);
   RadmonDetectorCarvedFlatVolumeConstructor &  operator=(const RadmonDetectorCarvedFlatVolumeConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorCarvedFlatVolumeConstructor.icc"
#endif /* RADMONDETECTORCARVEDFLATVOLUMECONSTRUCTOR_HH */
