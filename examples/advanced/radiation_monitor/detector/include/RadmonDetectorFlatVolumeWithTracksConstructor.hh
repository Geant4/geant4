//
// File name:     RadmonDetectorFlatVolumeWithTracksConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeWithTracksConstructor.hh,v 1.1 2005-09-27 13:55:51 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with tracks
//

#ifndef   RADMONDETECTORFLATVOLUMEWITHTRACKSCONSTRUCTOR_HH
 #define  RADMONDETECTORFLATVOLUMEWITHTRACKSCONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorVolumesWithTracksDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorFlatVolumeWithTracksConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorVolumesWithTracksDecorator<RadmonDetectorFlatVolumeComponent> >
 {
  public:
   inline                                       RadmonDetectorFlatVolumeWithTracksConstructor();
   inline virtual                              ~RadmonDetectorFlatVolumeWithTracksConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorFlatVolumeWithTracksConstructor(const RadmonDetectorFlatVolumeWithTracksConstructor & copy);
   RadmonDetectorFlatVolumeWithTracksConstructor & operator=(const RadmonDetectorFlatVolumeWithTracksConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorFlatVolumeWithTracksConstructor.icc"
#endif /* RADMONDETECTORFLATVOLUMEWITHTRACKSCONSTRUCTOR_HH */
