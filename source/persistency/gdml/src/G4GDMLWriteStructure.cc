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

G4ThreeVector G4GDMLWriteStructure::getRotation(const G4RotationMatrix& mat) {

   G4double x,y,z;

   G4double cosb = sqrt(mat.xx()*mat.xx()+mat.yx()*mat.yx());

   if (cosb > 16*FLT_EPSILON) {

      x = atan2(mat.zy(),mat.zz());
      y = atan2(-mat.zx(),cosb);
      z = atan2(mat.yx(),mat.xx());
   }
   else {

      x = atan2(-mat.yz(),mat.yy());
      y = atan2(-mat.zx(),cosb);
      z = 0.0;
   }

   x = -x/CLHEP::degree;
   y = -y/CLHEP::degree;
   z = -z/CLHEP::degree;

   return G4ThreeVector(x,y,z);
}

void G4GDMLWriteStructure::physvolWrite(xercesc::DOMElement* element,const G4VPhysicalVolume* const physvol) {

   xercesc::DOMElement* physvolElement = newElement("physvol");

   xercesc::DOMElement* volumerefElement = newElement("volumeref");
   volumerefElement->setAttributeNode(newAttribute("ref",physvol->GetLogicalVolume()->GetName()));
   physvolElement->appendChild(volumerefElement);

   G4ThreeVector pos = physvol->GetObjectTranslation();
   G4ThreeVector rot = getRotation(physvol->GetObjectRotationValue());

   if (pos.x() != 0.0 || pos.y() != 0.0 || pos.z() != 0.0) positionWrite(physvolElement,pos);
   if (rot.x() != 0.0 || rot.y() != 0.0 || rot.z() != 0.0) rotationWrite(physvolElement,rot);

   element->appendChild(physvolElement);
}

void G4GDMLWriteStructure::positionWrite(xercesc::DOMElement* element,const G4ThreeVector& pos) {

   xercesc::DOMElement* positionElement = newElement("position");
   positionElement->setAttributeNode(newAttribute("x",pos.x()));
   positionElement->setAttributeNode(newAttribute("y",pos.y()));
   positionElement->setAttributeNode(newAttribute("z",pos.z()));
   positionElement->setAttributeNode(newAttribute("unit","mm"));
   element->appendChild(positionElement);
}

void G4GDMLWriteStructure::rotationWrite(xercesc::DOMElement* element,const G4ThreeVector& rot) {

   xercesc::DOMElement* rotationElement = newElement("rotation");
   rotationElement->setAttributeNode(newAttribute("x",rot.x()));
   rotationElement->setAttributeNode(newAttribute("y",rot.y()));
   rotationElement->setAttributeNode(newAttribute("z",rot.z()));
   rotationElement->setAttributeNode(newAttribute("unit","deg"));
   element->appendChild(rotationElement);
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
