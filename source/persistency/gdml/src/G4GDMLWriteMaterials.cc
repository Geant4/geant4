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

void G4GDMLWriteMaterials::materialWrite(xercesc::DOMElement* element) {

   const G4MaterialTable* materialList = G4Material::GetMaterialTable();
   const G4int materialCount = materialList->size();

   for (G4int i=0;i<materialCount;i++) {
   
      const G4Material* materialPtr = (*materialList)[i];

      xercesc::DOMElement* materialElement = newElement("material");
    
      materialElement->setAttributeNode(newAttribute("name",materialPtr->GetName()));
      materialElement->setAttributeNode(newAttribute("Z",materialPtr->GetZ()));
    
      atomWrite(materialElement,materialPtr->GetA());    
    
      DWrite(materialElement,materialPtr->GetDensity());
      
      element->appendChild(materialElement);
   }
}

void G4GDMLWriteMaterials::materialsWrite(xercesc::DOMElement* element) {

   xercesc::DOMElement* materialsElement = newElement("materials");

   materialWrite(materialsElement);

   element->appendChild(materialsElement);
}
