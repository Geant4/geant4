//
// File name:     RadmonVDetectorLabelledEntityConstructor.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorLabelledEntityConstructor.cc,v 1.2 2005-09-21 14:56:43 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonVDetectorLabelledEntityConstructor.hh"
#include "RadmonMessenger.hh"
#include "RadmonTokenizer.hh"
#include "RadmonMaterialsManager.hh"
#include "G4UIcommand.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"
   
G4double                                        RadmonVDetectorLabelledEntityConstructor :: GetAttributeAsDouble(const G4String & attributeName, double defaultValue) const
{
 G4String str;
 
 str=GetAttribute(attributeName, "#");
 if (str=="#")
  return defaultValue;
  
 G4String args[1];
 if (!RadmonMessenger::ProcessArguments(str, 1, args))
  return defaultValue;
  
 return G4UIcommand::ConvertToDouble(args[0]);
}



G4double                                        RadmonVDetectorLabelledEntityConstructor :: GetAttributeAsMeasure(const G4String & attributeName, const char * category, double defaultValue) const
{
 G4String str;
 
 str=GetAttribute(attributeName, "#");
 if (str=="#")
  return defaultValue;
 
 G4String args[2];
 if (!RadmonMessenger::ProcessArguments(str, 2, args))
  return defaultValue;

 G4double unit(RadmonMessenger::GetUnit(args[1], category));
 if (unit<=0.)
  return defaultValue;
  
 return G4UIcommand::ConvertToDouble(args[0])*unit;
}



G4int                                           RadmonVDetectorLabelledEntityConstructor :: GetAttributeAsInteger(const G4String & attributeName, G4int defaultValue) const
{
 G4String str;
 
 str=GetAttribute(attributeName, "#");
 if (str=="#")
  return defaultValue;
  
 G4String args[1];
 if (!RadmonMessenger::ProcessArguments(str, 1, args))
  return defaultValue;
  
 return G4UIcommand::ConvertToInt(args[0]);
}





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

 G4String str;
 str=GetAttribute(attributeName, "#");
 if (str=="#")
 {
  G4String materialName(material->GetName());
  
  RadmonMaterialsManager * instance(RadmonMaterialsManager::Instance());
  
  visAttribute->SetColor(instance->GetMaterialColor(materialName));
  visAttribute->SetVisibility(instance->GetMaterialVisibility(materialName));
  if (instance->GetMaterialForceSolid(materialName))
   visAttribute->SetForceSolid(true);
  else if (instance->GetMaterialForceWireframe(materialName))
   visAttribute->SetForceWireframe(true);
 }
 else
 {
  G4String arg;
  RadmonTokenizer args(str);
  
  G4bool forceSolid(false);
  G4bool forceWireframe(false);
  G4bool visible(true);
  G4double alpha(1.);
  G4double red(1.);
  G4double green(1.);
  G4double blue(1.);
  
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
