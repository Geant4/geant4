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

void G4GDMLWriteStructure::physvolWrite(xercesc::DOMElement* volumeElement,const G4VPhysicalVolume* const physvol,const G4Transform3D& invR) {
/*
   G4Transform3D P(physvol->GetObjectRotationValue().inverse(),physvol->GetObjectTranslation());
   G4Transform3D R = volumeWrite(physvol->GetLogicalVolume());   
   G4Transform3D T = invR*P*R;

   HepGeom::Scale3D scale;
   HepGeom::Rotate3D rotate;
   HepGeom::Translate3D translate;

   T.getDecomposition(scale,rotate,translate);

   G4ThreeVector scl(scale(0,0),scale(1,1),scale(2,2));
   G4ThreeVector rot = getAngles(rotate.getRotation());
   G4ThreeVector pos = T.getTranslation();

   xercesc::DOMElement* physvolElement = newElement("physvol");
   volumeElement->appendChild(physvolElement);
   physvolElement->setAttributeNode(newAttribute("name",physvol->GetName()));

   xercesc::DOMElement* volumerefElement = newElement("volumeref");
   physvolElement->appendChild(volumerefElement);
   volumerefElement->setAttributeNode(newAttribute("ref",physvol->GetLogicalVolume()->GetName()));

   if (scl.x() != 1.0 || scl.y() != 1.0 || scl.z() != 1.0) scaleWrite(physvolElement,scl);
   if (rot.x() != 0.0 || rot.y() != 0.0 || rot.z() != 0.0) rotationWrite(physvolElement,rot);
   if (pos.x() != 0.0 || pos.y() != 0.0 || pos.z() != 0.0) positionWrite(physvolElement,pos);*/
}

G4Transform3D G4GDMLWriteStructure::volumeWrite(const G4LogicalVolume* volumePtr) {

   G4Transform3D R;

   volumeListNode* iter = last;
   
   while (iter != 0) {
   
      if (iter->volumePtr == volumePtr) { 
      
         if (iter->last == last) return R;
      
         iter->next->prev = iter->last->prev; // cut
         iter->last->prev->next = iter->next;
	 
	 iter->next = last; // paste at the end
	 last->prev = iter;
         iter->last->prev = 0;
         last = iter->last;
      
         return R;
      }
      
      iter = iter->next;
   }

//***************************************

   volumeListNode* newNode = new volumeListNode;
   newNode->volumePtr = volumePtr;
   newNode->volumeElement = newElement("volume");
   newNode->volumeElement->setAttributeNode(newAttribute("name",volumePtr->GetName()));

   newNode->next = last;
   newNode->prev = 0;
   
   if (last != 0) last->prev = newNode;
   
   last = newNode;
   
//***************************************

   G4cout << "Touching volume: " << volumePtr->GetName() << G4endl;

   G4Transform3D invR = R.inverse();

   const G4int daughterCount = volumePtr->GetNoDaughters();

   for (G4int i=0;i<daughterCount;i++) {
   
      const G4VPhysicalVolume* physvol = volumePtr->GetDaughter(i);
      const G4LogicalVolume* logvol = physvol->GetLogicalVolume();

      volumeWrite(logvol);
   }

   newNode->last = last;

   return R;   
}

void G4GDMLWriteStructure::structureWrite(xercesc::DOMElement* gdmlElement,const G4LogicalVolume* worldvol) {

   xercesc::DOMElement* structureElement = newElement("structure");
   gdmlElement->appendChild(structureElement);

   last = 0;

   volumeWrite(worldvol);

   G4cout << G4endl;

   volumeListNode* iter = last;
   
   while (iter != 0) {
   
      G4cout << "Printing volume: " << iter->volumePtr->GetName() << G4endl;
      
      structureElement->appendChild(iter->volumeElement);
      
      volumeListNode* temp = iter;
      iter = iter->next;
      delete temp;
   }
}
