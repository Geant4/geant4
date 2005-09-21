//
// File name:     RadmonDetectorFlatVolumeConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeConstructor.hh,v 1.1 2005-09-21 14:50:24 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box
//

#ifndef   RADMONDETECTORFLATVOLUMECONSTRUCTOR_HH
 #define  RADMONDETECTORFLATVOLUMECONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorFlatVolumeConstructor : public RadmonTDetectorLayerConstructor<RadmonDetectorFlatVolumeComponent>
 {
  public:
   inline                                       RadmonDetectorFlatVolumeConstructor();
   inline virtual                              ~RadmonDetectorFlatVolumeConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorFlatVolumeConstructor(const RadmonDetectorFlatVolumeConstructor & copy);
   RadmonDetectorFlatVolumeConstructor &         operator=(const RadmonDetectorFlatVolumeConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorFlatVolumeConstructor.icc"
#endif /* RADMONDETECTORFLATVOLUMECONSTRUCTOR_HH */
