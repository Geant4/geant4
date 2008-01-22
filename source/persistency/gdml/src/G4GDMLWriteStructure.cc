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
   element->appendChild(physvolElement);

   xercesc::DOMElement* volumerefElement = newElement("volumeref");
   volumerefElement->setAttributeNode(newAttribute("ref",physvol->GetLogicalVolume()->GetName()));
   physvolElement->appendChild(volumerefElement);

   G4Transform3D transform(physvol->GetObjectRotationValue(),physvol->GetObjectTranslation());

   G4VSolid* solidPtr = physvol->GetLogicalVolume()->GetSolid();

   if (const G4ReflectedSolid* reflectedPtr = dynamic_cast<const G4ReflectedSolid*>(solidPtr)) {   // Resolve reflected solid
      
      transform = transform*reflectedPtr->GetTransform3D();
      solidPtr = reflectedPtr->GetConstituentMovedSolid();
   }

   G4Scale3D scale;
   G4Rotate3D rotation;
   G4Translate3D translation;

   transform.getDecomposition(scale,rotation,translation);

   G4ThreeVector scl(scale.xx(),scale.yy(),scale.zz());
   G4ThreeVector rot = -getAngles(rotation.getRotation());
   G4ThreeVector pos(translation.dx(),translation.dy(),translation.dz());

   if (scl.x() != 1.0 || scl.y() != 1.0 || scl.z() != 1.0) scaleWrite(physvolElement,scl);
   if (rot.x() != 0.0 || rot.y() != 0.0 || rot.z() != 0.0) rotationWrite(physvolElement,rot);
   if (pos.x() != 0.0 || pos.y() != 0.0 || pos.z() != 0.0) positionWrite(physvolElement,pos);
}

void G4GDMLWriteStructure::volumeWrite(xercesc::DOMElement* element) {

   const G4LogicalVolumeStore* volumeList = G4LogicalVolumeStore::GetInstance();
   const size_t volumeCount = volumeList->size();

   for (size_t i=0;i<volumeCount;i++) {

      const G4LogicalVolume* volumePtr = (*volumeList)[i];

      xercesc::DOMElement* volumeElement = newElement("volume");
      element->appendChild(volumeElement);

      volumeElement->setAttributeNode(newAttribute("name",volumePtr->GetName()));

      xercesc::DOMElement* materialrefElement = newElement("materialref");
      materialrefElement->setAttributeNode(newAttribute("ref",volumePtr->GetMaterial()->GetName()));
      volumeElement->appendChild(materialrefElement);

      G4VSolid* solidPtr = volumePtr->GetSolid();

      if (const G4ReflectedSolid* reflectedPtr = dynamic_cast<const G4ReflectedSolid*>(solidPtr)) {   // Resolve reflected solid
      
         solidPtr = reflectedPtr->GetConstituentMovedSolid();
      }

      xercesc::DOMElement* solidrefElement = newElement("solidref");
      solidrefElement->setAttributeNode(newAttribute("ref",solidPtr->GetName()));
      volumeElement->appendChild(solidrefElement);

      const G4int daughterCount = volumePtr->GetNoDaughters();

      for (G4int j=0;j<daughterCount;j++) {
      
         physvolWrite(volumeElement,volumePtr->GetDaughter(j));
      }
   }
}

void G4GDMLWriteStructure::structureWrite(xercesc::DOMElement* element) {

   xercesc::DOMElement* structureElement = newElement("structure");
   element->appendChild(structureElement);

   volumeWrite(structureElement);
}
