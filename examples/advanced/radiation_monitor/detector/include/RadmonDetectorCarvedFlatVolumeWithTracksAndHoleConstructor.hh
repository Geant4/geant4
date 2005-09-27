//
// File name:     RadmonDetectorCarvedFlatVolumeWithTracksAndHoleConstructor.hh
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorCarvedFlatVolumeWithTracksAndHoleConstructor.hh,v 1.1 2005-09-27 13:55:51 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//
// Description:   Builds a box with carved borders
//

#ifndef   RADMONDETECTORCARVEDFLATVOLUMEWITHTRACKSANDHOLECONSTRUCTOR_HH
 #define  RADMONDETECTORCARVEDFLATVOLUMEWITHTRACKSANDHOLECONSTRUCTOR_HH

 // Include files
 #include "RadmonTDetectorLayerConstructor.hh"
 #include "RadmonTDetectorVolumesWithTracksDecorator.hh"
 #include "RadmonTDetectorVolumesWithHoleDecorator.hh"
 #include "RadmonTDetectorCarvedVolumesDecorator.hh"
 #include "RadmonDetectorFlatVolumeComponent.hh"
 
 class RadmonDetectorCarvedFlatVolumeWithTracksAndHoleConstructor : public RadmonTDetectorLayerConstructor<RadmonTDetectorVolumesWithTracksDecorator<RadmonTDetectorCarvedVolumesDecorator<RadmonTDetectorVolumesWithHoleDecorator<RadmonDetectorFlatVolumeComponent> > > >
 {
  public:
   inline                                       RadmonDetectorCarvedFlatVolumeWithTracksAndHoleConstructor();
   inline virtual                              ~RadmonDetectorCarvedFlatVolumeWithTracksAndHoleConstructor();

  private:
  // Hidden constructors and operators
                                                RadmonDetectorCarvedFlatVolumeWithTracksAndHoleConstructor(const RadmonDetectorCarvedFlatVolumeWithTracksAndHoleConstructor & copy);
   RadmonDetectorCarvedFlatVolumeWithTracksAndHoleConstructor & operator=(const RadmonDetectorCarvedFlatVolumeWithTracksAndHoleConstructor & copy);
 };
 
 // Inline implementations
 #include "RadmonDetectorCarvedFlatVolumeWithTracksAndHoleConstructor.icc"
#endif /* RADMONDETECTORCARVEDFLATVOLUMEWITHTRACKSANDHOLECONSTRUCTOR_HH */
