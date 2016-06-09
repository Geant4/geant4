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
// File name:     RadmonMaterialsManager.cc
// Creation date: Sep 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonMaterialsManager.cc,v 1.4 2006/06/29 16:16:51 gunter Exp $
// Tag:           $Name: geant4-09-01 $
//

// Include files
#include "RadmonMaterialsManager.hh"
#include "RadmonMaterialsMessenger.hh"
#include "RadmonMaterialsDumpStyle.hh"
#include "G4UnitsTable.hh"
#include "G4VisAttributes.hh"
#include "G4Color.hh"

#include <iomanip>





RadmonMaterialsManager *                        RadmonMaterialsManager :: Instance(void)
{
 if (!instance)
 {
  instance=new RadmonMaterialsManager;
  
  if (!instance)
   G4Exception("RadmonMaterialsManager::Instance: RadmonMaterialsManager singleton not allocated.");
 }
 
 return instance;
}





G4Element &                                     RadmonMaterialsManager :: CreateElement(const G4String & elementName, const G4String & symbol, G4double zEff, G4double aEff)
{
 G4Element * element(FindElement(elementName));
 
 if (element!=0)
 {
  G4cout << "RadmonMaterialsManager::CreateElement: Element \"" << elementName << "\" just exists." << G4endl;
  return (* element);
 }
 
 element=new G4Element(elementName, symbol, zEff, aEff);
 
 if (element==0)
 {
  G4String text("RadmonMaterialsManager::CreateElement: Element \"");
  text+=elementName;
  text+="\" not created.";
  
  G4Exception(text);
 }
 
 return (* element);
}



G4Element &                                     RadmonMaterialsManager :: GetElement(const G4String & elementName)
{
 G4Element * element(FindElement(elementName));
 
 if (element)
  return (* element);

 G4String text("RadmonMaterialsManager::GetElement: Element \"");
 text+=elementName;
 text+="\" not found.";

 G4Exception(text);
 return (* element);
}



G4bool                                          RadmonMaterialsManager :: ExistsElement(const G4String & elementName) const
{
 return const_cast<RadmonMaterialsManager *>(this)->FindElement(elementName)!=0;
}





void                                            RadmonMaterialsManager :: CreateMaterial(const G4String & materialName, G4double density, G4int nComponents)
{
 if (FindMaterial(materialName)!=0 || FindIncompleteMaterial(materialName)!=0)
 {
  G4cout << "RadmonMaterialsManager::CreateMaterial: Material \"" << materialName << "\" just exists." << G4endl;
  return;
 }

 if (nComponents==0)
 {
  G4String text("RadmonMaterialsManager::CreateMaterial: Material \"");
  text+=materialName;
  text+="\" must have at least one component.";
  
  G4Exception(text);
 }

 G4Material * material(new G4Material(materialName, density, nComponents));
 
 if (material==0)
 {
  G4String text("RadmonMaterialsManager::CreateMaterial: Material \"");
  text+=materialName;
  text+="\" not created.";
  
  G4Exception(text);
 }
 
 incompleteMaterialsList.push_back(material);
 G4VisAttributes * visAttributes(new G4VisAttributes);
 attributesMap[materialName]=visAttributes;
 visAttributes->SetColor(1., 1., 1., 1.);
 visAttributes->SetVisibility(true);
 
 return;
}



void                                            RadmonMaterialsManager :: AddComponentByAtoms(const G4String & materialName, const G4String & elementName, G4int nAtoms)
{
 G4Material * material(FindIncompleteMaterialOrAbort(materialName));
 G4Element & element(GetElement(elementName)); 
 
 material->AddElement(&element, nAtoms);
 
 UpdateIncompleteStatus(materialName);
}



void                                            RadmonMaterialsManager :: AddComponentByFraction(const G4String & materialName, const G4String & componentName, G4double fraction)
{
 G4Material * material(FindIncompleteMaterialOrAbort(materialName));
 
 G4Material * componentMat(FindMaterial(componentName));
 if (componentMat)
 {
  material->AddMaterial(componentMat, fraction);
  UpdateIncompleteStatus(materialName);
  return;
 }
 
 G4Element * componentElem(FindElement(componentName));
 if (componentElem)
 {
  material->AddElement(componentElem, fraction);
  UpdateIncompleteStatus(materialName);
  return;
 }
 
 G4String text("RadmonMaterialsManager::AddComponentByFraction: Component  \"");
 text+=componentName;
 text+="\" not found."; 
 
 G4Exception(text);
}



G4Material &                                    RadmonMaterialsManager :: GetMaterial(const G4String & materialName)
{
 G4Material *material(FindMaterial(materialName));
 
 if (material)
  return (* material);

 G4String text("RadmonMaterialsManager::GetMaterial: Material \"");
 text+=materialName;
 text+="\" not found.";

 G4Exception(text);
 return (* material);  
}



G4bool                                          RadmonMaterialsManager :: ExistsMaterial(const G4String & materialName) const
{
 return const_cast<RadmonMaterialsManager *>(this)->FindMaterial(materialName)!=0;
}
 
 



void                                            RadmonMaterialsManager :: SetMaterialColor(const G4String & materialName, const G4Color & color)
{
 if (attributesMap.find(materialName)==attributesMap.end())
 {
  G4cout << "RadmonMaterialsManager::SetMaterialColor: Material \"" << materialName << "\" does not exist." << G4endl;
  return;
 }

 attributesMap[materialName]->SetColor(color);
}



void                                            RadmonMaterialsManager :: SetMaterialVisibility(const G4String & materialName, G4bool visibility)
{
 if (attributesMap.find(materialName)==attributesMap.end())
 {
  G4cout << "RadmonMaterialsManager::SetMaterialVisibility: Material \"" << materialName << "\" does not exist." << G4endl;
  return;
 }

 attributesMap[materialName]->SetVisibility(visibility);
}



void                                            RadmonMaterialsManager :: SetMaterialForceWireframe(const G4String & materialName, G4bool force)
{
 if (attributesMap.find(materialName)==attributesMap.end())
 {
  G4cout << "RadmonMaterialsManager::SetMaterialForceWireframe: Material \"" << materialName << "\" does not exist." << G4endl;
  return;
 }

 attributesMap[materialName]->SetForceWireframe(force);
}



void                                            RadmonMaterialsManager :: SetMaterialForceSolid(const G4String & materialName, G4bool force)
{
 if (attributesMap.find(materialName)==attributesMap.end())
 {
  G4cout << "RadmonMaterialsManager::SetMaterialForceSolid: Material \"" << materialName << "\" does not exist." << G4endl;
  return;
 }

 attributesMap[materialName]->SetForceSolid(force);
}




 
const G4Color &                                 RadmonMaterialsManager :: GetMaterialColor(const G4String & materialName) const
{
 MaterialAttributes::const_iterator i(attributesMap.find(materialName));

 if (i==attributesMap.end())
 {
  static G4Color white(1., 1., 1., 1.);
  G4cout << "RadmonMaterialsManager::GetMaterialColor: Material \"" << materialName << "\" does not exist." << G4endl;

  return white;
 }

 return i->second->GetColor();
}



G4bool                                          RadmonMaterialsManager :: GetMaterialVisibility(const G4String & materialName) const
{
 MaterialAttributes::const_iterator i(attributesMap.find(materialName));

 if (i==attributesMap.end())
 {
  G4cout << "RadmonMaterialsManager::GetMaterialVisibility: Material \"" << materialName << "\" does not exist." << G4endl;
  return true;
 }

 return i->second->IsVisible();
}



G4bool                                          RadmonMaterialsManager :: GetMaterialForceWireframe(const G4String & materialName) const
{
 MaterialAttributes::const_iterator i(attributesMap.find(materialName));

 if (i==attributesMap.end())
 {
  G4cout << "RadmonMaterialsManager::GetMaterialForceWireframe: Material \"" << materialName << "\" does not exist." << G4endl;
  return false;
 }

 if (!i->second->IsForceDrawingStyle())
  return false;

 return i->second->GetForcedDrawingStyle()==G4VisAttributes::wireframe;
}



G4bool                                          RadmonMaterialsManager :: GetMaterialForceSolid(const G4String & materialName) const
{
 MaterialAttributes::const_iterator i(attributesMap.find(materialName));

 if (i==attributesMap.end())
 {
  G4cout << "RadmonMaterialsManager::GetMaterialForceSolid: Material \"" << materialName << "\" does not exist." << G4endl;
  return false;
 }

 if (!i->second->IsForceDrawingStyle())
  return false;

 return i->second->GetForcedDrawingStyle()==G4VisAttributes::solid;
}





void                                            RadmonMaterialsManager :: Dump(std::ostream & out, const G4String &indent) const
{
 G4int width(RADMONMATERIALDUMP_INDENT_WIDTH-indent.length());
 if (width<0)
  width=0;
  
 G4String indent2(indent);
 indent2.prepend("  "); 
 G4int width2(width-2);
 if (width2<0)
  width2=0;

 G4String indent3(indent2);
 indent3.prepend("  "); 
 G4int width3(width2-2);
 if (width3<0)
  width3=0;

 out << indent << "Elements:\n";
 
 G4int n(GetNElements());
 
 if (n==0)
  out << indent2 << "No elements defined.\n\n";
 else
 {
  const G4Element * element;
 
  for (G4int i(0); i<n; i++)
  {
   element=&GetElement(i);
   out << indent2 << "Element # " << i << '\n'
       << indent3 << std::setw(width3); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Name"; out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << element->GetName() << "\"\n"
       << indent3 << std::setw(width3); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Zeff"; out.setf(std::ostream::right, std::ostream::adjustfield); out << " = "   << std::setprecision(RADMONMATERIALDUMP_DOUBLE_PRECISION) << std::setw(RADMONMATERIALDUMP_DOUBLE_WIDTH) << element->GetZ() << '\n'
       << indent3 << std::setw(width3); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Aeff"; out.setf(std::ostream::right, std::ostream::adjustfield); out << " = "   << std::setprecision(RADMONMATERIALDUMP_DOUBLE_PRECISION) << std::setw(RADMONMATERIALDUMP_DOUBLE_WIDTH) << G4BestUnit(element->GetA(), "Molar Mass") << "\n\n";
  }
 } 

 out << indent << "Materials:\n";
 
 n=GetNMaterials();
 
 if (n==0)
  out << indent2 << "No materials defined.\n";
 else
 {
  const G4Material * material;
  G4int m;
 
  G4String indent4(indent3);
  indent4.prepend("  "); 
  G4int width4(width3-2);
  if (width4<0)
   width4=0;
  
  const G4double * massFractions;
  const G4double * nAtoms;

  for (G4int i(0); i<n; i++)
  {
   if (i>0)
    out << '\n';
  
   material=&GetMaterial(i);
   out << indent2 << "Material # " << i << '\n'
       << indent3 << std::setw(width3); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Name";    out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << material->GetName() << "\"\n"
       << indent3 << std::setw(width3); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Density"; out.setf(std::ostream::right, std::ostream::adjustfield); out << " = "   << std::setprecision(RADMONMATERIALDUMP_DOUBLE_PRECISION) << std::setw(RADMONMATERIALDUMP_DOUBLE_WIDTH) << G4BestUnit(material->GetDensity(), "Volumic Mass") << '\n';

   m=material->GetNumberOfElements();
   massFractions=material->GetFractionVector();
   nAtoms=material->GetVecNbOfAtomsPerVolume();
   
   for (G4int j(0); j<m; j++)
   {
    out <<                           indent4 << std::setw(width4); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Element";       out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << material->GetElement(j)->GetName() << "\"\n"
        << std::setw(indent4.length()) << "" << std::setw(width4); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Mass fraction"; out.setf(std::ostream::right, std::ostream::adjustfield); out << " = " << std::setprecision(RADMONMATERIALDUMP_PERCENT_PRECISION) << std::setw(RADMONMATERIALDUMP_PERCENT_WIDTH) << (massFractions[j]/perCent) << " %\n"
        << std::setw(indent4.length()) << "" << std::setw(width4); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Abundance";     out.setf(std::ostream::right, std::ostream::adjustfield); out << " = " << std::setprecision(RADMONMATERIALDUMP_PERCENT_PRECISION) << std::setw(RADMONMATERIALDUMP_PERCENT_WIDTH) << (nAtoms[j]/(perCent*material->GetTotNbOfAtomsPerVolume())) << " %\n";
   }
  }
 }

 if (!incompleteMaterialsList.empty())
 {
  out << '\n' << indent << "Materials to be completed:\n";
  
  MaterialsList::const_iterator i(incompleteMaterialsList.begin());
  MaterialsList::const_iterator end(incompleteMaterialsList.end());
   
  while (i!=end)
  {
   out << indent2 << std::setw(width2); out.setf(std::ostream::left, std::ostream::adjustfield); out << "Name";    out.setf(std::ostream::right, std::ostream::adjustfield); out << " = \"" << (*i)->GetName() << "\"\n";
   i++;
  }
 } 
}





G4bool                                          RadmonMaterialsManager :: Insert(std::istream & /* in */)
{
 // TO BE DONE
 G4cout << "RadmonMaterialsManager::Insert(): PLEASE CHECK" << G4endl;

 return false; 
}



G4bool                                          RadmonMaterialsManager :: Save(std::ostream & /* out */) const
{
 // TO BE DONE
 G4cout << "RadmonMaterialsManager::Save(): PLEASE CHECK" << G4endl;

 return false; 
}





inline G4Element *                              RadmonMaterialsManager :: FindElement(const G4String & elementName)
{
 G4int n(GetNElements());
 G4Element * element;
 
 while (n>0)
 {
  n--;
  element=&GetElement(n);
  
  if (element->GetName()==elementName)
   return element;
 }
 
 return 0;
}



inline G4Material *                             RadmonMaterialsManager :: FindMaterial(const G4String & materialName)
{
 G4int n(GetNMaterials());
 G4Material * material;
 
 while (n>0)
 {
  n--;
  material=&GetMaterial(n);
  
  if (material->GetName()==materialName)
   return material;
 }
 
 return 0;
}



G4Material *                                    RadmonMaterialsManager :: FindIncompleteMaterial(const G4String & materialName)
{
 MaterialsList::iterator i(incompleteMaterialsList.begin());
 MaterialsList::iterator end(incompleteMaterialsList.end());
 
 while (i!=end)
 {
  if ((*i)->GetName()==materialName)
   return (*i);
 
  i++;
 }
 
 return 0;
}



G4Material *                                    RadmonMaterialsManager :: FindIncompleteMaterialOrAbort(const G4String & materialName)
{
 G4Material * material(FindIncompleteMaterial(materialName));
 
 if (material)
  return material;
  
 if (&GetMaterial(materialName))
 {
  G4String text("RadmonMaterialsManager::FindIncompleteMaterialOrAbort: Material \"");
  text+=materialName;
  text+="\" does not need further components.";

  G4Exception(text);
 }
 
 return 0;
}



void                                            RadmonMaterialsManager :: UpdateIncompleteStatus(const G4String & materialName)
{
 if (!FindMaterial(materialName))
  return;

 MaterialsList::iterator i(incompleteMaterialsList.begin());
 MaterialsList::iterator end(incompleteMaterialsList.end());

 while (i!=end)
 {
  if ((*i)->GetName()==materialName)
  {
   incompleteMaterialsList.erase(i);
   return;
  }
  
  i++;
 }
}





inline                                          RadmonMaterialsManager :: RadmonMaterialsManager()
:
 messenger(new RadmonMaterialsMessenger(this))
{
 new G4UnitDefinition("g/mole",  "g/mole",  "Molar Mass", g/mole);
 new G4UnitDefinition("mg/mole", "mg/mole", "Molar Mass", mg/mole);
 new G4UnitDefinition("kg/mole", "kg/mole", "Molar Mass", kg/mole);
 
 new G4Material("RADMON_VACUUM", 1., 1.01*g/mole, universe_mean_density, kStateGas, 0.1*kelvin, 1.e-19*pascal);
}



                                                RadmonMaterialsManager :: ~RadmonMaterialsManager()
{
 if (messenger)
  delete messenger;
}





RadmonMaterialsManager *                        RadmonMaterialsManager :: instance(0);
