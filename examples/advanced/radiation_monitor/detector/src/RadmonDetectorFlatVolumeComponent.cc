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
// File name:     RadmonDetectorFlatVolumeComponent.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeComponent.cc,v 1.2.2.2.4.1 2009/08/11 14:20:35 gcosmo Exp $
// Tag:           $Name: geant4-09-02-patch-04 $
//

// Include files
#include "RadmonDetectorFlatVolumeComponent.hh"
#include "RadmonDetectorLayerVolumesList.hh"
#include "RadmonDetectorLayerVolumeItem.hh"
#include "RadmonVDetectorLabelledEntityConstructor.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
 
RadmonDetectorFlatVolumeComponent :: ~RadmonDetectorFlatVolumeComponent()
{
 delete visAttributes;
 delete box;
}

RadmonDetectorLayerVolumesList * RadmonDetectorFlatVolumeComponent :: GenerateVolumesList(void) 
{
 G4double width(owner->GetWidth());
 if (width<0)
  return 0;

 G4double height(owner->GetHeight());
 if (height<0)
  return 0;
  
 G4double thickness(owner->GetThickness());
 if (thickness<0)
  return 0;

 G4Material * material(owner->GetMaterial("Material"));  
 if (!material)
  return 0;
  
 visAttributes=owner->AllocateVisAttributes("VisAttributes", material);
 box=new G4Box("FlatVolume", width/2., height/2., thickness/2.);

 RadmonDetectorLayerVolumesList * list=new RadmonDetectorLayerVolumesList;
 RadmonDetectorLayerVolumeItem * item=list->AppendItem();

 item->SetSolid(box);
 item->SetAttributes(visAttributes);
 item->SetMaterial(material);
 item->SetSensitiveDetector(owner->AllocateSensitiveDetector("SensitiveDetector", ""));
 item->SetName("FlatVolume");
 
 return list;
}
