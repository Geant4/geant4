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
// File name:     RadmonDetectorLayout.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonDetectorLayout.cc,v 1.3.2.2 2006/06/29 16:13:59 gunter Exp $
// Tag:           $Name: geant4-09-02 $
//

// Include files
#include "RadmonDetectorLayout.hh"



                                                RadmonDetectorLayout :: RadmonDetectorLayout()
{
}



                                                RadmonDetectorLayout :: ~RadmonDetectorLayout()
{
}





void                                            RadmonDetectorLayout :: EnableEnvironment(void)
{
 if (environment.IsEnabled())
  return;

 environment.Enable();
 NotifyChange();
}



void                                            RadmonDetectorLayout :: DisableEnvironment(void)
{
 if (!environment.IsEnabled())
  return;

 environment.Disable();
 NotifyChange();
}



G4bool                                          RadmonDetectorLayout :: IsEnabledEnvironment(void) const
{
 return environment.IsEnabled();
}





void                                            RadmonDetectorLayout :: SetEnvironmentType(const G4String & type)
{
 if (environment.GetType()==type)
  return;
 
 environment.SetType(type);
 
 if (environment.IsEnabled())
  NotifyChange();
}



const G4String &                                RadmonDetectorLayout :: GetEnvironmentType() const
{
 return environment.GetType();
}



G4int                                           RadmonDetectorLayout :: GetEnvironmentNAttributes(void) const
{
 return environment.GetNAttributes();
}



const G4String &                                RadmonDetectorLayout :: GetEnvironmentAttributeName(G4int index) const
{
 return environment.GetAttributeName(index);
}



void                                            RadmonDetectorLayout :: SetEnvironmentAttribute(const G4String & attributeName, const G4String & attributeValue)
{
 if (attributeName=="")
 {
  G4cout << "RadmonDetectorLayout::FindPlacement: \"\" is not a valid attribute name." << G4endl;
  return;
 }

 if (environment.GetAttribute(attributeName, attributeValue+"#")==attributeValue)
  return;
 
 environment.SetAttribute(attributeName, attributeValue);

 if (environment.IsEnabled())
  NotifyChange();
}



const G4String                                  RadmonDetectorLayout :: GetEnvironmentAttribute(const G4String & attributeName, const G4String & defaultAttributeValue) const
{
 return environment.GetAttribute(attributeName, defaultAttributeValue);
}



void                                            RadmonDetectorLayout :: ClearEnvironmentAttribute(const G4String & attributeName)
{
 if (!environment.ExistsAttribute(attributeName))
 {
  G4cout << "RadmonDetectorLayout::ClearEnvironmentAttribute: Attribute \"" << attributeName << "\" missing in environment." << G4endl;
  return;
 }
  
 environment.ClearAttribute(attributeName);

 if (environment.IsEnabled())
  NotifyChange();
}





void                                            RadmonDetectorLayout :: CreateMultilayer(const G4String & multilayerLabel)
{
 if (multilayersCollection.ExistsMultilayerByLabel(multilayerLabel))
 {
  G4cout << "RadmonDetectorLayout::CreateMultilayer: Multilayer \"" << multilayerLabel << "\" just exists." << G4endl;
  return;
 }
  
 RadmonDetectorMultilayerLayout & multilayer(multilayersCollection.CreateMultilayer());
 multilayer.SetLabel(multilayerLabel);
}



void                                            RadmonDetectorLayout :: RemoveMultilayer(const G4String & multilayerLabel)
{
 if (IsPlaced(multilayerLabel))
 {
  G4cout << "RadmonDetectorLayout::RemoveMultilayer: Multilayer \"" << multilayerLabel << "\" is placed and cannot be deleted. Remove the placement first." << G4endl;
  return;
 }

 if (!multilayersCollection.ExistsMultilayerByLabel(multilayerLabel))
 {
  G4cout << "RadmonDetectorLayout::RemoveMultilayer: Multilayer \"" << multilayerLabel << "\" does not exist." << G4endl;
  return;
 }
 
 multilayersCollection.RemoveMultilayersByLabel(multilayerLabel);
}



void                                            RadmonDetectorLayout :: SetMultilayerWidth(const G4String & multilayerLabel, G4double width)
{
 RadmonDetectorMultilayerLayout * multilayer(FindMultilayer(multilayerLabel));
 
 if (!multilayer)
  return;
 
 if (multilayer->GetWidth()==width)
  return;
 
 multilayer->SetWidth(width);
 
 if (IsPlaced(multilayerLabel))
  NotifyChange();
}



G4double                                        RadmonDetectorLayout :: GetMultilayerWidth(const G4String & multilayerLabel) const
{
 const RadmonDetectorMultilayerLayout * multilayer(FindMultilayer(multilayerLabel));
 
 if (!multilayer)
  return 0.;
  
 return multilayer->GetWidth();
}



void                                            RadmonDetectorLayout :: SetMultilayerHeight(const G4String & multilayerLabel, G4double height)
{
 RadmonDetectorMultilayerLayout * multilayer(FindMultilayer(multilayerLabel));
 
 if (!multilayer)
  return;
 
 if (multilayer->GetHeight()==height)
  return;
 
 multilayer->SetHeight(height);
 
 if (IsPlaced(multilayerLabel))
  NotifyChange();
}



G4double                                        RadmonDetectorLayout :: GetMultilayerHeight(const G4String & multilayerLabel) const
{
 const RadmonDetectorMultilayerLayout * multilayer(FindMultilayer(multilayerLabel));

 if (!multilayer)
  return 0.;

 return multilayer->GetHeight();
}



G4double                                        RadmonDetectorLayout :: GetMultilayerTotalThickness(const G4String & multilayerLabel) const
{
 const RadmonDetectorMultilayerLayout * multilayer(FindMultilayer(multilayerLabel));

 if (!multilayer)
  return 0.;

 return multilayer->GetTotalThickness();
}





void                                            RadmonDetectorLayout :: AppendLayerToMultilayer(const G4String & multilayerLabel, const G4String & layerLabel)
{
 RadmonDetectorMultilayerLayout * multilayer(FindMultilayer(multilayerLabel));
 
 if (!multilayer)
  return;
 
 if (multilayer->ExistsLayerByLabel(layerLabel))
 {
  G4cout << "RadmonDetectorLayout::AppendLayerToMultilayer: Layer \"" << layerLabel << "\" just exists in multilayer \"" << multilayerLabel << "\"." << G4endl;
  return;
 }
  
 RadmonDetectorLayerLayout & layer(multilayer->AppendLayer());
 
 layer.SetLabel(layerLabel);
 
 if (IsPlaced(multilayerLabel))
  NotifyChange();
}



void                                            RadmonDetectorLayout :: RemoveLayerFromMultilayer(const G4String & multilayerLabel, const G4String & layerLabel)
{
 RadmonDetectorMultilayerLayout * multilayer(FindMultilayer(multilayerLabel));
 
 if (!multilayer)
  return;
 
 if (!multilayer->ExistsLayerByLabel(layerLabel))
 {
  G4cout << "RadmonDetectorLayout::RemoveLayerFromMultilayer: Layer \"" << layerLabel << "\" does not exist in multilayer \"" << multilayerLabel << "\"." << G4endl;
  return;
 }
  
 multilayer->RemoveLayersByLabel(layerLabel);

 if (IsPlaced(multilayerLabel))
  NotifyChange();
}



void                                            RadmonDetectorLayout :: RemoveAllLayersFromMultilayer(const G4String & multilayerLabel)
{
 RadmonDetectorMultilayerLayout * multilayer(FindMultilayer(multilayerLabel));
 
 if (!multilayer)
  return;
  
 multilayer->RemoveAllLayers();

 if (IsPlaced(multilayerLabel))
  NotifyChange();
}



G4int                                           RadmonDetectorLayout :: GetMultilayerNLayers(const G4String & multilayerLabel) const
{
 const RadmonDetectorMultilayerLayout * multilayer(FindMultilayer(multilayerLabel));
 
 if (!multilayer)
  return 0;
  
 return multilayer->GetNLayers();
}



const G4String &                                RadmonDetectorLayout :: GetMultilayerLayerLabel(const G4String & multilayerLabel, G4int index) const
{
 const RadmonDetectorMultilayerLayout * multilayer(FindMultilayer(multilayerLabel));
 
 if (!multilayer)
  return GetNullStr();
  
 return multilayer->GetLayer(index).GetLabel();
}





void                                            RadmonDetectorLayout :: SetLayerThickness(const G4String & multilayerLabel, const G4String & layerLabel, G4double thickness)
{
 RadmonDetectorLayerLayout * layer(FindLayer(multilayerLabel, layerLabel));
 
 if (!layer)
  return;
  
 if (layer->GetThickness()==thickness)
  return;
 
 layer->SetThickness(thickness);
 
 if (IsPlaced(multilayerLabel))
  NotifyChange();
}



G4double                                        RadmonDetectorLayout :: GetLayerThickness(const G4String & multilayerLabel, const G4String & layerLabel) const
{
 const RadmonDetectorLayerLayout * layer(FindLayer(multilayerLabel, layerLabel));
 
 if (!layer)
  return 0.;
 
 return layer->GetThickness();
}



void                                            RadmonDetectorLayout :: SetLayerType(const G4String & multilayerLabel, const G4String & layerLabel, const G4String & type)
{
 RadmonDetectorLayerLayout * layer(FindLayer(multilayerLabel, layerLabel));
 
 if (!layer)
  return;
  
 if (layer->GetType()==type)
  return;
 
 layer->SetType(type);
 
 if (IsPlaced(multilayerLabel))
  NotifyChange();
}



const G4String &                                RadmonDetectorLayout :: GetLayerType(const G4String & multilayerLabel, const G4String & layerLabel) const
{
 const RadmonDetectorLayerLayout * layer(FindLayer(multilayerLabel, layerLabel));
 
 if (!layer)
  return GetNullStr();
  
 return layer->GetType();
}





G4int                                           RadmonDetectorLayout :: GetLayerNAttributes(const G4String & multilayerLabel, const G4String & layerLabel) const
{
 const RadmonDetectorLayerLayout * layer(FindLayer(multilayerLabel, layerLabel));
 
 if (!layer)
  return 0;

 return layer->GetNAttributes();
}



const G4String &                                RadmonDetectorLayout :: GetLayerAttributeName(const G4String & multilayerLabel, const G4String & layerLabel, G4int index) const
{
 const RadmonDetectorLayerLayout * layer(FindLayer(multilayerLabel, layerLabel));
 
 if (!layer)
  return GetNullStr();

 return layer->GetAttributeName(index);
}



void                                            RadmonDetectorLayout :: SetLayerAttribute(const G4String & multilayerLabel, const G4String & layerLabel, const G4String & attributeName, const G4String & attributeValue)
{
 if (attributeName=="")
 {
  G4cout << "RadmonDetectorLayout::FindPlacement: \"\" is not a valid attribute name." << G4endl;
  return;
 }

 RadmonDetectorLayerLayout * layer(FindLayer(multilayerLabel, layerLabel));
 
 if (!layer)
  return;
  
 if (layer->GetAttribute(attributeName, attributeValue+"#")==attributeValue)
  return;
 
 layer->SetAttribute(attributeName, attributeValue);
 
 if (IsPlaced(multilayerLabel))
  NotifyChange();
}



const G4String                                  RadmonDetectorLayout :: GetLayerAttribute(const G4String & multilayerLabel, const G4String & layerLabel, const G4String & attributeName, const G4String & defaultAttributeValue) const
{
 const RadmonDetectorLayerLayout * layer(FindLayer(multilayerLabel, layerLabel));
 
 if (!layer)
  return G4String("");
  
 return layer->GetAttribute(attributeName, defaultAttributeValue);
}



void                                            RadmonDetectorLayout :: ClearLayerAttribute(const G4String & multilayerLabel, const G4String & layerLabel, const G4String & attributeName)
{
 RadmonDetectorLayerLayout * layer(FindLayer(multilayerLabel, layerLabel));
 
 if (!layer)
  return;
  
 if (!layer->ExistsAttribute(attributeName))
 {
  G4cout << "RadmonDetectorLayout::ClearLayerAttribute: Attribute \"" << attributeName << "\" missing in layer \"" << layerLabel << "\" of multilayer \"" << multilayerLabel << "\"." << G4endl;
  return;
 }
  
 layer->ClearAttribute(attributeName);
 
 if (IsPlaced(multilayerLabel))
  NotifyChange();
}





void                                            RadmonDetectorLayout :: CreatePlacement(const G4String & placementLabel, const G4String & multilayerName)
{
 if (multilayerPlacementsCollection.ExistsPlacementByLabel(placementLabel))
 {
  G4cout << "RadmonDetectorLayout::CreateMultilayer: Placement \"" << placementLabel << "\" just exists." << G4endl;
  return;
 }
  
 RadmonDetectorMultilayerPlacementLayout & placement(multilayerPlacementsCollection.CreatePlacement());
 placement.SetLabel(placementLabel);
 placement.SetMultilayerLabel(multilayerName);

 NotifyChange();
}



G4int                                           RadmonDetectorLayout :: GetNPlacements() const
{
 return multilayerPlacementsCollection.GetNPlacements();
}



const G4String &                                RadmonDetectorLayout :: GetPlacementLabel(G4int index) const
{
 return multilayerPlacementsCollection.GetPlacement(index).GetLabel();
}



void                                            RadmonDetectorLayout :: RemovePlacement(const G4String & placementLabel)
{
 if (!multilayerPlacementsCollection.ExistsPlacementByLabel(placementLabel))
 {
  G4cout << "RadmonDetectorLayout::CreateMultilayer: Placement \"" << placementLabel << "\" does not exist." << G4endl;
  return;
 }
  
 multilayerPlacementsCollection.RemovePlacementByLabel(placementLabel);
 NotifyChange();
}





const G4String &                                RadmonDetectorLayout :: GetPlacementMultilayerType(const G4String & placementLabel) const
{
 const RadmonDetectorMultilayerPlacementLayout * placement(FindPlacement(placementLabel));

 if (!placement)
  return GetNullStr();
  
 return placement->GetMultilayerLabel();
}



void                                            RadmonDetectorLayout :: SetPlacementPosition(const G4String & placementLabel, const G4ThreeVector & position)
{
 RadmonDetectorMultilayerPlacementLayout * placement(FindPlacement(placementLabel));

 if (!placement)
  return;

 placement->SetAbsolutePosition(position);
 NotifyChange();  
}



const G4ThreeVector &                           RadmonDetectorLayout :: GetPlacementPosition(const G4String & placementLabel) const
{
 const RadmonDetectorMultilayerPlacementLayout * placement(FindPlacement(placementLabel));

 if (!placement)
  return GetNullPosition();
  
 return placement->GetAbsolutePosition();
}



void                                            RadmonDetectorLayout :: SetPlacementPosition(const G4String & placementLabel, const G4String & originLabel, const G4ThreeVector & offset)
{
 RadmonDetectorMultilayerPlacementLayout * placement(FindPlacement(placementLabel));

 if (!placement)
  return; 
  
 const RadmonDetectorMultilayerPlacementLayout * origin(FindPlacement(originLabel));
 if (!origin)
  return;
  
 placement->SetRelativePosition(*origin, offset);
 NotifyChange();   
}



void                                            RadmonDetectorLayout :: SetPlacementRotation(const G4String & placementLabel, const G4RotationMatrix & rotation)
{
 RadmonDetectorMultilayerPlacementLayout * placement(FindPlacement(placementLabel));

 if (!placement)
  return;

 placement->SetAbsoluteRotation(rotation);
 NotifyChange();  
}



const G4RotationMatrix &                        RadmonDetectorLayout :: GetPlacementRotation(const G4String & placementLabel) const
{
 const RadmonDetectorMultilayerPlacementLayout * placement(FindPlacement(placementLabel));

 if (!placement)
  return GetNullRotationMatrix();
  
 return placement->GetAbsoluteRotation();
}



void                                            RadmonDetectorLayout :: SetPlacementRotation(const G4String & placementLabel, const G4String & originLabel, const G4RotationMatrix & relativeRotation)
{
 RadmonDetectorMultilayerPlacementLayout * placement(FindPlacement(placementLabel));

 if (!placement)
  return; 
  
 const RadmonDetectorMultilayerPlacementLayout * origin(FindPlacement(originLabel));
 if (!origin)
  return;
  
 placement->SetRelativeRotation(*origin, relativeRotation);
 NotifyChange();   
}





void                                            RadmonDetectorLayout :: DumpLayout(std::ostream & out) const
{
 const G4String indent("  - ");

 out << "- Environment\n";
 environment.DumpLayout(out, indent);
 
 out << "\n- Multilayers\n";
 multilayersCollection.DumpLayout(out, indent);
 
 out << "\n- Placements\n";
 multilayerPlacementsCollection.DumpLayout(out, indent);
}





G4bool                                          RadmonDetectorLayout :: Load(std::istream & /*in*/)
{
 // TO BE DONE
 G4cout << "RadmonDetectorLayout::Load(): PLEASE CHECK" << G4endl;

 return false; 
}



G4bool                                          RadmonDetectorLayout :: Save(std::ostream & /*out*/) const
{
 // TO BE DONE
 G4cout << "RadmonDetectorLayout::Save(): PLEASE CHECK" << G4endl;

 return false; 
}





inline bool                                     RadmonDetectorLayout :: IsPlaced(const G4String & multilayerLabel)
{
 G4int n(multilayerPlacementsCollection.GetNPlacements());
 
 while (n>0)
 {
  n--;
  
  if (multilayerPlacementsCollection.GetPlacement(n).GetMultilayerLabel()==multilayerLabel)
   return true;
 }
 
 return false;
}



inline RadmonDetectorMultilayerLayout *         RadmonDetectorLayout :: FindMultilayer(const G4String & multilayerLabel)
{
 if (!multilayersCollection.ExistsMultilayerByLabel(multilayerLabel))
 {
  G4cout << "RadmonDetectorLayout::FindMultilayer: Multilayer \"" << multilayerLabel << "\" does not exist." << G4endl;
  return 0;
 }

 return &multilayersCollection.FindMultilayerByLabel(multilayerLabel);
}



inline const RadmonDetectorMultilayerLayout *   RadmonDetectorLayout :: FindMultilayer(const G4String & multilayerLabel) const
{
 return const_cast<RadmonDetectorLayout *>(this)->FindMultilayer(multilayerLabel);
}



inline RadmonDetectorMultilayerPlacementLayout * RadmonDetectorLayout :: FindPlacement(const G4String & placementLabel)
{
 if (!multilayerPlacementsCollection.ExistsPlacementByLabel(placementLabel))
 {
  G4cout << "RadmonDetectorLayout::FindPlacement: Placement \"" << placementLabel << "\" does not exist." << G4endl;
  return 0;
 }

 return & multilayerPlacementsCollection.FindPlacementByLabel(placementLabel);
}



inline const RadmonDetectorMultilayerPlacementLayout * RadmonDetectorLayout :: FindPlacement(const G4String & placementLabel) const
{
 return const_cast<RadmonDetectorLayout *>(this)->FindPlacement(placementLabel);
}



inline RadmonDetectorLayerLayout *              RadmonDetectorLayout :: FindLayer(const G4String & multilayerLabel, const G4String & layerLabel)
{
 RadmonDetectorMultilayerLayout * multilayer(FindMultilayer(multilayerLabel));
 
 if (!multilayer)
  return 0;
 
 if (!multilayer->ExistsLayerByLabel(layerLabel))
 {
  G4cout << "RadmonDetectorLayout::FindLayer: Layer \"" << layerLabel << "\" does not exist in multilayer \"" << multilayerLabel << "\"." << G4endl;
  return 0;
 }
 
 return &multilayer->FindLayerByLabel(layerLabel);
}



inline const RadmonDetectorLayerLayout *        RadmonDetectorLayout :: FindLayer(const G4String & multilayerLabel, const G4String & layerLabel) const
{
 return const_cast<RadmonDetectorLayout *>(this)->FindLayer(multilayerLabel, layerLabel);
}





inline G4String &                               RadmonDetectorLayout :: GetNullStr() const
{
 static G4String *nullStr(0);
 
 if (nullStr==0)
  nullStr=new G4String("");
  
 return *nullStr;
}



inline G4ThreeVector &                          RadmonDetectorLayout :: GetNullPosition() const
{
 static G4ThreeVector *nullPosition(0);
 
 if (nullPosition==0)
  nullPosition=new G4ThreeVector();
  
 return *nullPosition;
}



inline G4RotationMatrix &                       RadmonDetectorLayout :: GetNullRotationMatrix() const
{
 static G4RotationMatrix * nullRotationMatrix(0);
 
 if (nullRotationMatrix==0)
  nullRotationMatrix=new G4RotationMatrix();
  
 return *nullRotationMatrix;
}
