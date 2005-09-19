//
// File name:     RadmonVDetectorLabelledEntityConstructor.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonVDetectorLabelledEntityConstructor.cc,v 1.1 2005-09-19 19:38:41 capra Exp $
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
   
G4double                                        RadmonVDetectorLabelledEntityConstructor :: GetAttributeAsDouble(const G4String & attributeName, double defaultValue)
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



G4double                                        RadmonVDetectorLabelledEntityConstructor :: GetAttributeAsMeasure(const G4String & attributeName, const char * category, double defaultValue)
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



G4int                                           RadmonVDetectorLabelledEntityConstructor :: GetAttributeAsInteger(const G4String & attributeName, G4int defaultValue)
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





G4VisAttributes *                               RadmonVDetectorLabelledEntityConstructor :: AllocateVisAttributes(const G4String & attributeName, const G4String & materialName)
{
 G4VisAttributes * visAttribute=new G4VisAttributes;

 G4String str;
 str=GetAttribute(attributeName, "#");
 if (str=="#")
 {
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
  G4double red(1.);
  G4double green(1.);
  G4double blue(1.);
  
  G4bool forceFound(false);
  G4bool visibleFound(false);
  G4int missingColors(3);
  
  for (G4int i(0); i<5; i++)
  {
   if (args.eos())
    break;
    
   arg=args();
   
   if (arg=="hidden" || arg=="visible")
   {
    if (visibleFound || missingColors==1 || missingColors==2)
    {
     missingColors=-1;
     break;
    }
    
    visibleFound=true;
    
    visible=(arg=="visible");
   }
   else if (arg=="wireframe" || arg=="solid" || arg=="default")
   {
    if (forceFound || missingColors==1 || missingColors==2)
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
      case 3:
       red=component;
       break;

      case 2:
       green=component;
       break;

      case 1:
       blue=component;
       break;
     }

     missingColors--;
    }
   }
  }
  
  if (missingColors!=0 && missingColors!=3)
   G4cout << "RadmonVDetectorLabelledEntityConstructor::AllocateVisAttribute: Invalid visualization attribute in \"" << attributeName << "\"." << G4endl; 
  else
  {
   visAttribute->SetColor(red, green, blue);
   visAttribute->SetVisibility(visible);
   if (forceWireframe)
    visAttribute->SetForceWireframe(true);
   else if (forceSolid)
    visAttribute->SetForceSolid(true);
  }
 }
 
 return visAttribute;
}
