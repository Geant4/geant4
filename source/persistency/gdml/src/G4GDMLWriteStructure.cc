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

#include "G4GDMLWriteStructure.hh"

void G4GDMLWriteStructure::physvolWrite(xercesc::DOMElement* element,const G4VPhysicalVolume* const physvol) {

   xercesc::DOMElement* physvolElement = newElement("physvol");

   xercesc::DOMElement* volumerefElement = newElement("volumeref");
   volumerefElement->setAttributeNode(newAttribute("ref",physvol->GetLogicalVolume()->GetName()));

   physvolElement->appendChild(volumerefElement);
   element->appendChild(physvolElement);
}

void G4GDMLWriteStructure::volumeWrite(xercesc::DOMElement* element) {

   const G4LogicalVolumeStore* volumeList = G4LogicalVolumeStore::GetInstance();
   const G4int volumeCount = volumeList->size();

   for (G4int i=0;i<volumeCount;i++) {

      const G4LogicalVolume* volumePtr = (*volumeList)[i];

      xercesc::DOMElement* volumeElement = newElement("volume");
      volumeElement->setAttributeNode(newAttribute("name",volumePtr->GetName()));

      xercesc::DOMElement* materialrefElement = newElement("materialref");
      materialrefElement->setAttributeNode(newAttribute("ref",volumePtr->GetMaterial()->GetName()));
      volumeElement->appendChild(materialrefElement);

      xercesc::DOMElement* solidrefElement = newElement("solidref");
      solidrefElement->setAttributeNode(newAttribute("ref",volumePtr->GetSolid()->GetName()));
      volumeElement->appendChild(solidrefElement);

      const G4int daughterCount = volumePtr->GetNoDaughters();

      for (G4int j=0;j<daughterCount;j++) {
      
         physvolWrite(volumeElement,volumePtr->GetDaughter(j));
      }

      element->appendChild(volumeElement);
   }
}

void G4GDMLWriteStructure::structureWrite(xercesc::DOMElement* element) {

   xercesc::DOMElement* structureElement = newElement("structure");

   volumeWrite(structureElement);

   element->appendChild(structureElement);
}
