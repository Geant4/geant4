//
// File name:     RadmonDetectorSimpleBoxConstructor.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorSimpleBoxConstructor.cc,v 1.1 2005-09-19 19:39:52 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorSimpleBoxConstructor.hh"
#include "RadmonMaterialsManager.hh"
#include "RadmonMessenger.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4UIcommand.hh"
#include "globals.hh"
 
                                                RadmonDetectorSimpleBoxConstructor :: ~RadmonDetectorSimpleBoxConstructor()
{
 delete logicalVolume;
 delete visAttributes;
 delete box;
}


   
G4LogicalVolume *                               RadmonDetectorSimpleBoxConstructor :: ConstructLogicalVolume(void)
{
 if (logicalVolume)
  G4Exception("RadmonDetectorSimpleBoxConstructor::ConstructLogicalVolume: Called twice.");

 G4double width(GetAttributeAsMeasure("_WIDTH", "Length", -1.));
 if (width<0)
 {
  width=GetAttributeAsMeasure("Width", "Length", -1.);
  
  if (width<0)
  {
   G4cout << "RadmonDetectorSimpleBoxConstructor::ConstructLogicalVolume: \"Width\" attribute not defined." << G4endl;
   return 0;
  }
 }

 G4double height(GetAttributeAsMeasure("_HEIGHT", "Length", -1.));
 if (height<0)
 {
  height=GetAttributeAsMeasure("Height", "Length", -1.);
  
  if (height<0)
  {
   G4cout << "RadmonDetectorSimpleBoxConstructor::ConstructLogicalVolume: \"Height\" attribute not defined." << G4endl;
   return 0;
  }
 }
 
 G4double thickness(GetAttributeAsMeasure("_THICKNESS", "Length", -1.));
 if (thickness<0)
 {
  thickness=GetAttributeAsMeasure("Thickness", "Length", -1.);
  
  if (thickness<0)
  {
   G4cout << "RadmonDetectorSimpleBoxConstructor::ConstructLogicalVolume: \"Thickness\" attribute not defined." << G4endl;
   return 0;
  }
 }
  
 G4String materialStr(GetAttribute("Material", "#")); 
 if (materialStr=="#")
 {
  G4cout << "RadmonDetectorSimpleBoxConstructor::ConstructLogicalVolume: \"Material\" attribute not defined." << G4endl;
  return 0;
 }
 
 RadmonMaterialsManager * manager(RadmonMaterialsManager::Instance());
 
 if (!manager->ExistsMaterial(materialStr))
 {
  G4cout << "RadmonDetectorSimpleBoxConstructor::ConstructLogicalVolume: \"" << materialStr << "\" material not existent." << G4endl;
  return 0;
 }
 
 G4Material & material(manager->GetMaterial(materialStr));
 visAttributes=AllocateVisAttributes("VisAttributes", materialStr);

 box=new G4Box(GetLabel(), width/2., thickness/2., height/2.);
 logicalVolume=new G4LogicalVolume(box, &material, GetLabel()+"LV", 0, 0, 0);
 logicalVolume->SetVisAttributes(visAttributes);
 
 return logicalVolume;
}



RadmonVDetectorLabelledEntityConstructor *      RadmonDetectorSimpleBoxConstructor :: New(void)
{
 return new RadmonDetectorSimpleBoxConstructor;
}
