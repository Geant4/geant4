//
// File name:     RadmonDetectorCarvedFlatVolumeWithTracksConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorCarvedFlatVolumeWithTracksConstructor.hh,v 1.1 2005-09-27 13:55:51 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with carved borders and tracks
//

#ifndef   RADMONDETECTORCARVEDFLATVOLUMEWITHTRACKSCONSTRUCTOR_HH
 #define  RADMONDETECTORCARVEDFLATVOLUMEWITHTRACKSCONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorCarvedVolumesDecorator.hh"
 #include "RadmonTDetectorVolumesWithTracksDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorCarvedFlatVolumeWithTracksConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorVolumesWithTracksDecorator<RadmonTDetectorCarvedVolumesDecorator<RadmonDetectorFlatVolumeComponent> > >
 {
  public:
   inline                                       RadmonDetectorCarvedFlatVolumeWithTracksConstructor();
   inline virtual                              ~RadmonDetectorCarvedFlatVolumeWithTracksConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorCarvedFlatVolumeWithTracksConstructor(const RadmonDetectorCarvedFlatVolumeWithTracksConstructor & copy);
   RadmonDetectorCarvedFlatVolumeWithTracksConstructor & operator=(const RadmonDetectorCarvedFlatVolumeWithTracksConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorCarvedFlatVolumeWithTracksConstructor.icc"
#endif /* RADMONDETECTORCARVEDFLATVOLUMEWITHTRACKSCONSTRUCTOR_HH */
