//
// File name:     RadmonDetectorConstruction.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorConstruction.cc,v 1.3 2005-09-19 19:42:13 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonDetectorConstruction.hh"

#include "RadmonVDetectorLayout.hh"
#include "RadmonVDetectorEntityConstructor.hh"
#include "RadmonVDetectorEntitiesConstructorsFactory.hh"
#include "RadmonMaterialsManager.hh"

#include "G4RunManager.hh"
#include "G4UIcommand.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"



                                                RadmonDetectorConstruction :: RadmonDetectorConstruction(RadmonVDetectorLayout * layout, RadmonVDetectorEntitiesConstructorsFactory * factory)
:
 detectorLayout(layout),
 constructorsFactory(factory),
 environmentConstructor(0),
 environmentPhysicalVolume(0),
 environmentLogicalVolume(0),
 environmentSolid(0)
{
 if (detectorLayout==0)
  G4Exception("RadmonDetectorConstruction::RadmonDetectorConstruction: layout==0.");

 if (factory==0)
  G4Exception("RadmonDetectorConstruction::RadmonDetectorConstruction: factory==0.");
 
 detectorLayout->AttachObserver(this);
}



                                                RadmonDetectorConstruction :: ~RadmonDetectorConstruction()
{
 detectorLayout->DetachObserver(this);

 Destruct();
 
 delete constructorsFactory;
}





G4VPhysicalVolume *                             RadmonDetectorConstruction :: Construct(void)
{
 // This method will not change the detector layout
 RadmonVDetectorLayout const * const & const_detectorLayout(detectorLayout);

 if (const_detectorLayout->IsEnabledEnvironment())
  BuildEnvironmentFromType(const_detectorLayout->GetEnvironmentType());
 else
  BuildEnvironmentSphere();

 environmentPhysicalVolume=new G4PVPlacement(0, G4ThreeVector(), "worldPV", environmentLogicalVolume, 0, false, 0);
 
 G4int placement(const_detectorLayout->GetNPlacements());
 while (placement>0)
 {
  placement--;
  BuildMultilayer(placement);
 }
 
 return environmentPhysicalVolume;
}





void                                            RadmonDetectorConstruction :: OnLayoutChange(void)
{
 Destruct();
 G4RunManager::GetRunManager()->DefineWorldVolume(Construct(), true);
}




void                                            RadmonDetectorConstruction :: Destruct(void)
{
 // Layers destruction
 LayerItem item;
 
 while (!layersStack.empty())
 {
  item=layersStack.top();
  
  environmentLogicalVolume->RemoveDaughter(item.second);
  delete item.second;
  delete item.first;
  
  layersStack.pop();
 }

 // Environment destruction
 delete environmentPhysicalVolume;
 environmentPhysicalVolume=0;

 if (environmentConstructor)
 {
  delete environmentConstructor;
  environmentConstructor=0;
 }
 else
 {
  delete environmentLogicalVolume;
  delete environmentSolid;
  environmentSolid=0;
 }

 environmentLogicalVolume=0;
}





void                                            RadmonDetectorConstruction :: BuildEnvironmentFromType(const G4String & type)
{
 environmentConstructor=constructorsFactory->CreateEntityConstructor(type);
 
 if (environmentConstructor==0)
 {
  G4cout << "RadmonDetectorConstruction::BuildEnvironmentFromType: Environment type \"" << type << "\" not found. Environment ignored." << G4endl;
  BuildEnvironmentSphere();
 }

 // This method will not change the detector layout
 RadmonVDetectorLayout const * const & const_detectorLayout(detectorLayout);
 
 G4int n(const_detectorLayout->GetEnvironmentNAttributes());
 G4String name;
 G4String value;
 
 while (n>0)
 {
  n--;
  name=const_detectorLayout->GetEnvironmentAttributeName(n);
  value=const_detectorLayout->GetEnvironmentAttribute(name, G4String());
  
  environmentConstructor->SetEntityAttribute(name, value);
 }
 
 environmentLogicalVolume=environmentConstructor->ConstructLogicalVolume();
 
 if (environmentLogicalVolume==0)
 {
  delete environmentConstructor;
  environmentConstructor=0;
  G4cout << "RadmonDetectorConstruction::BuildEnvironmentFromType: Environment type \"" << type << "\" not built. Environment ignored." << G4endl;
  BuildEnvironmentSphere();
 }
}



void                                            RadmonDetectorConstruction :: BuildEnvironmentSphere(void)
{
 // This method will not change the detector layout
 RadmonVDetectorLayout const * const & const_detectorLayout(detectorLayout);
 
 G4double maxRadius(-0.1);
 G4double radius;
 G4ThreeVector vect;
 
 G4int placement(const_detectorLayout->GetNPlacements());
 G4String placementLabel;
 G4String multilayerLabel;
 
 while (placement>0)
 {
  placement--;
  placementLabel=const_detectorLayout->GetPlacementLabel(placement);
  multilayerLabel=const_detectorLayout->GetPlacementMultilayerType(placementLabel);
  
  vect.setX(const_detectorLayout->GetMultilayerWidth(multilayerLabel));
  vect.setY(const_detectorLayout->GetMultilayerHeight(multilayerLabel));
  vect.setZ(const_detectorLayout->GetMultilayerTotalThickness(multilayerLabel));
  radius=const_detectorLayout->GetPlacementPosition(placementLabel).getR()+vect.getR()/2.;
 
  if (radius>maxRadius)
   maxRadius=radius;
 }
 
 if (maxRadius<=0.)
  maxRadius=5*m;

 G4Material & vacuum(RadmonMaterialsManager::Instance()->GetMaterial("RADMON_VACUUM"));
 
 environmentSolid=new G4Sphere("worldSphere", 0.*mm, maxRadius, 0.*deg, 360.*deg, 0.*deg, 180.*deg);
 environmentLogicalVolume=new G4LogicalVolume(environmentSolid, &vacuum, "worldLV", 0, 0, 0);
}



void                                            RadmonDetectorConstruction :: BuildMultilayer(G4int index)
{
 // This method will not change the detector layout
 RadmonVDetectorLayout const * const & const_detectorLayout(detectorLayout);
 
 G4String placementLabel(const_detectorLayout->GetPlacementLabel(index));
 G4String multilayerLabel(const_detectorLayout->GetPlacementMultilayerType(placementLabel));
 G4RotationMatrix rotation(const_detectorLayout->GetPlacementRotation(placementLabel));
 G4ThreeVector position(const_detectorLayout->GetPlacementPosition(placementLabel));
  
 G4double width(const_detectorLayout->GetMultilayerWidth(multilayerLabel)); 
 if (width<=0.)
  return;

 G4String widthStr(G4UIcommand::ConvertToString(width/mm)+" mm");

 G4double height(const_detectorLayout->GetMultilayerHeight(multilayerLabel)); 
 if (height<=0.)
  return;

 G4String heightStr(G4UIcommand::ConvertToString(height/mm)+" mm");
 
 G4double y(const_detectorLayout->GetMultilayerTotalThickness(multilayerLabel)); 
 if (y<=0.)
  return;
  
 y/=2.;
 
 G4int n(const_detectorLayout->GetMultilayerNLayers(multilayerLabel));
 if (n==0)
  return;
  
 while (n>0)
 {
  n--;
  
  G4String layerLabel(const_detectorLayout->GetMultilayerLayerLabel(multilayerLabel, n));
  G4String layerType(const_detectorLayout->GetLayerType(multilayerLabel, layerLabel));

  G4double layerThickness(const_detectorLayout->GetLayerThickness(multilayerLabel, layerLabel));
  if (layerThickness<=0.)
   continue;
   
  RadmonVDetectorEntityConstructor * entityConstructor(constructorsFactory->CreateEntityConstructor(layerType));
  
  if (entityConstructor==0)
  {
   G4cout << "RadmonDetectorConstruction::BuildMultilayer: Layer type \"" << layerType << "\" not found. Layer ignored." << G4endl;
   continue;
  }

  G4int m(const_detectorLayout->GetLayerNAttributes(multilayerLabel, layerLabel));
  G4String name;
  G4String value;
 
  while (m>0)
  {
   m--;
   name=const_detectorLayout->GetLayerAttributeName(multilayerLabel, layerLabel, m);
   value=const_detectorLayout->GetLayerAttribute(multilayerLabel, layerLabel, name, G4String());
  
   entityConstructor->SetEntityAttribute(name, value);
  }
  
  entityConstructor->SetEntityAttribute("_WIDTH", widthStr);
  entityConstructor->SetEntityAttribute("_HEIGHT", heightStr);
  entityConstructor->SetEntityAttribute("_THICKNESS", G4UIcommand::ConvertToString(layerThickness/mm)+" mm");
  
  layerThickness/=2.;
  
  G4LogicalVolume * entityLogicalVolume(entityConstructor->ConstructLogicalVolume());
  
  if (entityLogicalVolume==0)
  {
   delete entityConstructor;
   G4cout << "RadmonDetectorConstruction::BuildMultilayer: Layer \"" << layerLabel << "\" from \"" << multilayerLabel << "\" not built. Layer ignored." << G4endl;
   continue;
  }
  
  y-=layerThickness;
  G4ThreeVector layerPosition(rotation.colY());
  layerPosition*=y;
  layerPosition+=position;
  
  G4VPhysicalVolume * entityPhysicalVolume(new G4PVPlacement(G4Transform3D(rotation, layerPosition), entityLogicalVolume, placementLabel+"_"+multilayerLabel+"_"+layerLabel, environmentLogicalVolume, false, 0));
  
  y-=layerThickness;
  
  layersStack.push(LayerItem(entityConstructor, entityPhysicalVolume));
 } 
}
