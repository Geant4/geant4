//
// File name:     RadmonDetectorLayerVolumeItem.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumeItem.cc,v 1.2 2005-11-25 01:53:30 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorLayerVolumeItem.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Material.hh"
#include "G4SDManager.hh"

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

 if (volumeSensitiveDetector)
 {
  G4SDManager * manager(G4SDManager::GetSDMpointer());

  if (manager)
   if (! manager->FindSensitiveDetector(volumeSensitiveDetector->GetFullPathName(), false))
    delete volumeSensitiveDetector;
 }
}





G4LogicalVolume *                               RadmonDetectorLayerVolumeItem :: GetLogicalVolume(void)
{
 if (! volumeLogical)
 {
  if (volumeSensitiveDetector)
  {
   G4SDManager * manager(G4SDManager::GetSDMpointer());
   
   if (manager)
   {
    G4VSensitiveDetector * sensitiveDetector(manager->FindSensitiveDetector(volumeSensitiveDetector->GetFullPathName(), false));
  
    if (!sensitiveDetector)
     manager->AddNewDetector(volumeSensitiveDetector);
    else if (sensitiveDetector!=volumeSensitiveDetector)
    {
     delete volumeSensitiveDetector;
     volumeSensitiveDetector=sensitiveDetector;
    }
  
    volumeSensitiveDetector->Activate(true);
   }
  }
  
  volumeLogical=new G4LogicalVolume(volumeSolid, volumeMaterial, volumeName+"LV", 0, volumeSensitiveDetector, 0);
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
