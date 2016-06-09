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
// File name:     RadmonVDetectorLabelledEntityConstructor.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorLabelledEntityConstructor.cc,v 1.7 2006/06/29 16:14:13 gunter Exp $
// Tag:           $Name: geant4-08-02 $
//

// Include files
#include "RadmonVDetectorLabelledEntityConstructor.hh"
#include "RadmonTokenizer.hh"
#include "RadmonMaterialsManager.hh"
#include "RadmonSensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4UIcommand.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
   
G4double                                        RadmonVDetectorLabelledEntityConstructor :: GetWidth(void) const
{
 G4double value(GetAttributeAsMeasure("_WIDTH", "Length", -1.));
 
 if (value>=0)
  return value;

 value=GetAttributeAsMeasure("Width", "Length", -1.);
 
 if (value>=0)
  return value;

 G4cout << "RadmonVDetectorLabelledEntityConstructor::GetWidth: \"Width\" attribute not defined." << G4endl;
 return -1;
}



G4double                                        RadmonVDetectorLabelledEntityConstructor :: GetHeight(void) const
{
 G4double value(GetAttributeAsMeasure("_HEIGHT", "Length", -1.));
 
 if (value>=0)
  return value;

 value=GetAttributeAsMeasure("Height", "Length", -1.);
 
 if (value>=0)
  return value;

 G4cout << "RadmonVDetectorLabelledEntityConstructor::GetHeight: \"Height\" attribute not defined." << G4endl;
 return -1;
}



G4double                                        RadmonVDetectorLabelledEntityConstructor :: GetThickness(void) const
{
 G4double value(GetAttributeAsMeasure("_THICKNESS", "Length", -1.));
 
 if (value>=0)
  return value;

 value=GetAttributeAsMeasure("Thickness", "Length", -1.);
 
 if (value>=0)
  return value;

 G4cout << "RadmonVDetectorLabelledEntityConstructor::GetThickness: \"Thickness\" attribute not defined." << G4endl;
 return -1;
}





G4Material *                                    RadmonVDetectorLabelledEntityConstructor :: GetMaterial(const G4String & attributeName) const
{
 G4String materialStr(GetAttribute(attributeName, "#")); 
 if (materialStr=="#")
 {
  G4cout << "RadmonDetectorSimpleBoxConstructor::ConstructLogicalVolume: \"" << attributeName << "\" attribute not defined." << G4endl;
  return 0;
 }
 
 RadmonMaterialsManager * manager(RadmonMaterialsManager::Instance()); 
 if (!manager->ExistsMaterial(materialStr))
 {
  G4cout << "RadmonDetectorSimpleBoxConstructor::ConstructLogicalVolume: \"" << materialStr << "\" material not existent." << G4endl;
  return 0;
 }
 
 return &manager->GetMaterial(materialStr);
}





G4VisAttributes *                               RadmonVDetectorLabelledEntityConstructor :: AllocateVisAttributes(const G4String & attributeName, const G4Material * material) const
{
 G4VisAttributes * visAttribute=new G4VisAttributes;

 G4String materialName(material->GetName());
  
 RadmonMaterialsManager * instance(RadmonMaterialsManager::Instance());
  
 visAttribute->SetColor(instance->GetMaterialColor(materialName));
 visAttribute->SetVisibility(instance->GetMaterialVisibility(materialName));
 if (instance->GetMaterialForceSolid(materialName))
  visAttribute->SetForceSolid(true);
 else if (instance->GetMaterialForceWireframe(materialName))
  visAttribute->SetForceWireframe(true);
  
 G4String str;
 str=GetAttribute(attributeName, "#");
 if (str!="#")
 {
  G4String arg;
  RadmonTokenizer args(str);
  
  G4bool forceSolid(visAttribute->GetForcedDrawingStyle()==G4VisAttributes::solid && visAttribute->IsForceDrawingStyle());
  G4bool forceWireframe(visAttribute->GetForcedDrawingStyle()==G4VisAttributes::wireframe && visAttribute->IsForceDrawingStyle());
  G4bool visible(visAttribute->IsVisible());
  G4double alpha(visAttribute->GetColor().GetAlpha());
  G4double red(visAttribute->GetColor().GetRed());
  G4double green(visAttribute->GetColor().GetGreen());
  G4double blue(visAttribute->GetColor().GetBlue());
  
  G4bool forceFound(false);
  G4bool visibleFound(false);
  G4int missingColors(4);
  
  for (G4int i(0); i<6; i++)
  {
   if (args.eos())
    break;
    
   arg=args();
   
   if (arg=="hidden" || arg=="visible")
   {
    if (visibleFound || missingColors==2 || missingColors==3)
    {
     missingColors=-1;
     break;
    }
    
    visibleFound=true;
    
    visible=(arg=="visible");
   }
   else if (arg=="wireframe" || arg=="solid" || arg=="default")
   {
    if (forceFound || missingColors==2 || missingColors==3)
    {
     missingColors=-1;
     break;
    }
    
    forceFound=true;
    
    forceWireframe=(arg=="wireframe");
    forceSolid=(arg=="solid");
   }
   else
   {
    G4double component(G4UIcommand::ConvertToDouble(arg));
    if (component>1 || component<0 || missingColors==0)
    {
     missingColors=-1;
     break;
    }
    else
    {
     switch(missingColors)
     {
      case 4:
       red=component;
       break;

      case 3:
       green=component;
       break;

      case 2:
       blue=component;
       break;

      case 1:
       alpha=component;
       break;
     }

     missingColors--;
    }
   }
  }
  
  if (missingColors!=0 && missingColors!=4 && missingColors!=1)
   G4cout << "RadmonVDetectorLabelledEntityConstructor::AllocateVisAttribute: Invalid visualization attribute in \"" << attributeName << "\"." << G4endl; 
  else
  {
   visAttribute->SetColor(red, green, blue, alpha);
   visAttribute->SetVisibility(visible);
   if (forceWireframe)
    visAttribute->SetForceWireframe(true);
   else if (forceSolid)
    visAttribute->SetForceSolid(true);
  }
 }
 
 return visAttribute;
}





G4VSensitiveDetector *                          RadmonVDetectorLabelledEntityConstructor :: AllocateSensitiveDetector(const G4String & attributeName, const G4String & defaultAttrbuteValue) const
{
 G4String name(GetAttribute(attributeName, defaultAttrbuteValue));
 
 if (name.empty())
  return 0;
  
 G4SDManager * manager(G4SDManager::GetSDMpointer());
  
 G4VSensitiveDetector * sensitiveDetector(manager->FindSensitiveDetector(name, false));

 if (sensitiveDetector)
  return sensitiveDetector;
  
 return new RadmonSensitiveDetector(name);
}

