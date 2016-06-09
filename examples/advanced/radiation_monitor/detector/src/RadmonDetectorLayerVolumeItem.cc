//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// File name:     RadmonDetectorLayerVolumeItem.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayerVolumeItem.cc,v 1.2.2.2.4.1 2009/08/11 14:20:35 gcosmo Exp $
// Tag:           $Name: geant4-09-02-patch-04 $
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

G4LogicalVolume * RadmonDetectorLayerVolumeItem :: GetLogicalVolume(void)
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
