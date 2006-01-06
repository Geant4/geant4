//
// File name:     RadmonApplicationDetectorSetup.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationDetectorSetup.cc,v 1.2 2006-01-06 12:52:31 guatelli Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonApplicationDetectorSetup.hh"
#include "RadmonApplicationOptions.hh"

#include "RadmonDetectorLabelledEntitiesConstructorsFactory.hh"

#include "RadmonDetectorFlatVolumeConstructor.hh"
#include "RadmonDetectorFlatVolumeWithHoleConstructor.hh"
#include "RadmonDetectorFlatVolumeWithGroundConstructor.hh"
#include "RadmonDetectorFlatVolumeWithGroundAndKeyMarksConstructor.hh"
#include "RadmonDetectorFlatVolumeWithGroundAndKeyMarksAndHoleConstructor.hh"
#include "RadmonDetectorFlatVolumeWithTracksConstructor.hh"
#include "RadmonDetectorFlatVolumeWithPinsConstructor.hh"
#include "RadmonDetectorFlatVolumeWithPadsConstructor.hh"
#include "RadmonDetectorFlatVolumeWithTracksAndHoleConstructor.hh"
#include "RadmonDetectorCarvedFlatVolumeWithHoleConstructor.hh"
#include "RadmonDetectorCarvedFlatVolumeConstructor.hh"
#include "RadmonDetectorCarvedFlatVolumeWithGroundConstructor.hh"
#include "RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor.hh"
#include "RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksAndHoleConstructor.hh"
#include "RadmonDetectorCarvedFlatVolumeWithTracksConstructor.hh"
#include "RadmonDetectorCarvedFlatVolumeWithPinsConstructor.hh"
#include "RadmonDetectorCarvedFlatVolumeWithTracksAndHoleConstructor.hh"



#define DECLARE_DETECTOR_CONSTRUCTOR(name)      constructor=new name();                                                                  \
                                                if (constructor==0)                                                                      \
                                                {                                                                                        \
                                                 G4cerr << currentOptions.ApplicationName() << ": Cannot allocate " #name "." << G4endl; \
                                                 return false;                                                                           \
                                                }                                                                                        \
                                                factory->AppendLabelledEntityConstructor(constructor)

G4bool RadmonApplicationDetectorSetup :: CreateDetectorEntityConstructors(RadmonDetectorLabelledEntitiesConstructorsFactory * factory)
{
 RadmonVDetectorLabelledEntityConstructor * constructor;
 
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeWithHoleConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeWithGroundConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeWithGroundAndKeyMarksConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeWithGroundAndKeyMarksAndHoleConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeWithTracksConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeWithPinsConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeWithTracksAndHoleConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorFlatVolumeWithPadsConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorCarvedFlatVolumeConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorCarvedFlatVolumeWithHoleConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorCarvedFlatVolumeWithTracksConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorCarvedFlatVolumeWithPinsConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorCarvedFlatVolumeWithTracksAndHoleConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorCarvedFlatVolumeWithGroundConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksConstructor);
 DECLARE_DETECTOR_CONSTRUCTOR(RadmonDetectorCarvedFlatVolumeWithGroundAndKeyMarksAndHoleConstructor);
 
 return true;
}
