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

#include "G4GDMLWriteSolids.hh"

void G4GDMLWriteSolids::boxWrite(xercesc::DOMElement* solidsElement,const G4Box* const box) {

   xercesc::DOMElement* boxElement = newElement("box");
   boxElement->setAttributeNode(newAttribute("name",box->GetName()));
   boxElement->setAttributeNode(newAttribute("x",2.0*box->GetXHalfLength()));
   boxElement->setAttributeNode(newAttribute("y",2.0*box->GetYHalfLength()));
   boxElement->setAttributeNode(newAttribute("z",2.0*box->GetZHalfLength()));
   boxElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(boxElement);
}

void G4GDMLWriteSolids::tessellatedWrite(xercesc::DOMElement* solidsElement,const G4TessellatedSolid* const tessellated) {

   xercesc::DOMElement* tessellatedElement = newElement("tessellated");
   tessellatedElement->setAttributeNode(newAttribute("name",tessellated->GetName()));

   const size_t n = tessellated->GetNumberOfFacets();
   
   for (size_t i = 0;i<n;i++) {
   
      const G4VFacet* facet = tessellated->GetFacet(i);
   }

   solidsElement->appendChild(tessellatedElement);
}

void G4GDMLWriteSolids::solidsWrite(xercesc::DOMElement* element) {

   const G4SolidStore* solidList = G4SolidStore::GetInstance();
   const G4int solidCount = solidList->size();

   xercesc::DOMElement* solidsElement = newElement("solids");

   for (G4int i=0;i<solidCount;i++) {
   
      const G4VSolid* solidPtr = (*solidList)[i];

      if (const G4Box* boxPtr = dynamic_cast<const G4Box*>(solidPtr)) { boxWrite(solidsElement,boxPtr); } else
      if (const G4TessellatedSolid* tessellatedPtr = dynamic_cast<const G4TessellatedSolid*>(solidPtr)) { tessellatedWrite(solidsElement,tessellatedPtr); }
   }

   element->appendChild(solidsElement);
}
