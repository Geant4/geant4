//
// File name:     RadmonDetectorLayerVolumeItem.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumeItem.cc,v 1.1 2005-09-21 14:52:32 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorLayerVolumeItem.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"

                                                RadmonDetectorLayerVolumeItem :: ~RadmonDetectorLayerVolumeItem()
{
 if (volumePhysical)
 {
  volumeMother->GetLogicalVolume()->RemoveDaughter(volumePhysical);
  delete volumePhysical;
 }
 
 if (volumeLogical->GetNoDaughters())
  G4Exception("RadmonDetectorLayerVolumeItem::~RadmonDetectorLayerVolumeItem: Logical volume removed prior of daughters");
 
 delete volumeLogical;
}





G4LogicalVolume *                               RadmonDetectorLayerVolumeItem :: GetLogicalVolume(void)
{
 if (! volumeLogical)
 {
  volumeLogical=new G4LogicalVolume(volumeSolid, volumeMaterial, volumeName+"LV", 0, 0, 0);
  if (volumeAttributes)
    volumeLogical->SetVisAttributes(volumeAttributes);
  
  if (volumeMother)
   volumePhysical=new G4PVPlacement(&volumeRotation, volumePosition, volumeLogical, volumeName+"LV", volumeMother->GetLogicalVolume(), false, 0);
 }

 return volumeLogical; 
}



void                                            RadmonDetectorLayerVolumeItem :: Assertion(void)
{
 if (volumeLogical)
  G4Exception("RadmonDetectorLayerVolumeItem::Assertion: You cannot change a RadmonDetectorLayerVolumeItem after GetLogicalVolume call.");
}
