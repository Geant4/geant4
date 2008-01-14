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

void G4GDMLWriteSolids::boxWrite(xercesc::DOMElement* element,const G4Box* const box) {

   xercesc::XMLString::transcode("box",tempStr,99);
   xercesc::DOMElement* boxElement = doc->createElement(tempStr);

   xercesc::XMLString::transcode("name",tempStr,99);
   xercesc::DOMAttr* name = doc->createAttribute(tempStr);
   
   xercesc::XMLString::transcode("zoli",tempStr,99);
   name->setValue(tempStr);
   
   boxElement->setAttributeNode(name);
   
   element->appendChild(boxElement);
}

void G4GDMLWriteSolids::solidsWrite(xercesc::DOMElement* element) {

   const G4SolidStore* solidList = G4SolidStore::GetInstance();
   const G4int solidCount = solidList->size();

   xercesc::XMLString::transcode("solids",tempStr,99);
   xercesc::DOMElement* solids = doc->createElement(tempStr);
   element->appendChild(solids);

   for (G4int i=0;i<solidCount;i++) {
   
      const G4VSolid* solidPtr = (*solidList)[i];

      if (const G4Box* boxPtr = dynamic_cast<const G4Box*>(solidPtr)) { boxWrite(solids,boxPtr); }
   }
}
