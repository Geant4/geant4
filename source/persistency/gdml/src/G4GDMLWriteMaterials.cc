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
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLWriteMaterials.hh"

void G4GDMLWriteMaterials::atomWrite(xercesc::DOMElement* element,G4double a) {

   a /= (g/mole);

   xercesc::DOMElement* atomElement = newElement("atom");
   atomElement->setAttributeNode(newAttribute("unit","g/mole"));
   atomElement->setAttributeNode(newAttribute("value",a));
   element->appendChild(atomElement);
}

void G4GDMLWriteMaterials::DWrite(xercesc::DOMElement* element,G4double d) {

   d /= (g/cm3);

   xercesc::DOMElement* DElement = newElement("D");
   DElement->setAttributeNode(newAttribute("unit","g/cm3"));
   DElement->setAttributeNode(newAttribute("value",d));
   element->appendChild(DElement);
}

void G4GDMLWriteMaterials::PWrite(xercesc::DOMElement* element,G4double P) {

   P /= pascal;

   xercesc::DOMElement* PElement = newElement("P");
   PElement->setAttributeNode(newAttribute("unit","pascal"));
   PElement->setAttributeNode(newAttribute("value",P));
   element->appendChild(PElement);
}

void G4GDMLWriteMaterials::TWrite(xercesc::DOMElement* element,G4double T) {

   T /= kelvin;

   xercesc::DOMElement* TElement = newElement("T");
   TElement->setAttributeNode(newAttribute("unit","K"));
   TElement->setAttributeNode(newAttribute("value",T));
   element->appendChild(TElement);
}

void G4GDMLWriteMaterials::isotopeWrite(xercesc::DOMElement* element,const G4Isotope* const isotopePtr) {

   xercesc::DOMElement* isotopeElement = newElement("isotope");
   element->appendChild(isotopeElement);

   isotopeElement->setAttributeNode(newAttribute("name",isotopePtr->GetName()));
}

void G4GDMLWriteMaterials::elementWrite(xercesc::DOMElement* element,const G4Element* const elementPtr) {

   const G4MaterialTable* materialList = G4Material::GetMaterialTable();
   const size_t materialCount = materialList->size();

   // If the element is registered both as a material and an element, we dump it only once as a material!

   for (size_t i=0;i<materialCount;i++)
      if ((*materialList)[i]->GetName() == elementPtr->GetName()) return;

   xercesc::DOMElement* elementElement = newElement("element");
   element->appendChild(elementElement);

   elementElement->setAttributeNode(newAttribute("name",elementPtr->GetName()));

   elementElement->setAttributeNode(newAttribute("Z",elementPtr->GetZ()));
   atomWrite(elementElement,elementPtr->GetA());
}

void G4GDMLWriteMaterials::materialWrite(xercesc::DOMElement* element,const G4Material* const materialPtr) {

   xercesc::DOMElement* materialElement = newElement("material");
   element->appendChild(materialElement);

   G4State state = materialPtr->GetState();

   G4String state_str("undefined");
   
   if (state==kStateSolid) { state_str = "solid"; } else
   if (state==kStateLiquid) { state_str = "liquid"; } else
   if (state==kStateGas) { state_str = "gas"; }

   materialElement->setAttributeNode(newAttribute("name",materialPtr->GetName()));
   materialElement->setAttributeNode(newAttribute("state",state_str));

   DWrite(materialElement,materialPtr->GetDensity());
   PWrite(materialElement,materialPtr->GetPressure());
   TWrite(materialElement,materialPtr->GetTemperature());
  
   const size_t NumberOfElements = materialPtr->GetNumberOfElements();

   if (NumberOfElements>1) { 
   
      const G4double* MassFractionVector = materialPtr->GetFractionVector();

      for (size_t i=0;i<NumberOfElements;i++) {
      
         xercesc::DOMElement* fractionElement = newElement("fraction");
         materialElement->appendChild(fractionElement);

         fractionElement->setAttributeNode(newAttribute("n",MassFractionVector[i]));
         fractionElement->setAttributeNode(newAttribute("ref",materialPtr->GetElement(i)->GetName()));
      }
 
   } else {
   
      materialElement->setAttributeNode(newAttribute("Z",materialPtr->GetZ()));
      atomWrite(materialElement,materialPtr->GetA());
   }
}

void G4GDMLWriteMaterials::materialsWrite(xercesc::DOMElement* element) {

   G4cout << "Writing materials..." << G4endl;

   xercesc::DOMElement* materialsElement = newElement("materials");
   element->appendChild(materialsElement);

   const G4IsotopeTable* isotopeList = G4Isotope::GetIsotopeTable();
   const size_t isotopeCount = isotopeList->size();

   for (size_t i=0;i<isotopeCount;i++)
      isotopeWrite(materialsElement,(*isotopeList)[i]);

   const G4ElementTable* elementList = G4Element::GetElementTable();
   const size_t elementCount = elementList->size();

   for (size_t i=0;i<elementCount;i++)
      elementWrite(materialsElement,(*elementList)[i]);
   
   const G4MaterialTable* materialList = G4Material::GetMaterialTable();
   const size_t materialCount = materialList->size();

   for (size_t i=0;i<materialCount;i++)
      materialWrite(materialsElement,(*materialList)[i]);
}
