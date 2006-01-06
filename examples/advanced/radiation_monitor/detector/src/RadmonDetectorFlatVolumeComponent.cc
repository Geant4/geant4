//
// File name:     RadmonDetectorFlatVolumeComponent.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorFlatVolumeComponent.cc,v 1.3 2006-01-06 12:52:32 guatelli Exp $
// Tag:           $Name: not supported by cvs2svn $
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
