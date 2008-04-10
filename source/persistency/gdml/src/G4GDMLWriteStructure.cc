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
   if (pos.x() != 0.0 || pos.y() != 0.0 || pos.z() != 0.0) positionWrite(physvolElement,pos);
}

G4Transform3D G4GDMLWriteStructure::volumeWrite(const G4LogicalVolume* volumePtr) {

   G4Transform3D R;
   G4VSolid* solidPtr = volumePtr->GetSolid();

   if (G4ReflectedSolid* refl = dynamic_cast<G4ReflectedSolid*>(solidPtr)) {
   
      R = refl->GetTransform3D();
      solidPtr = refl->GetConstituentMovedSolid();
   }

   for (int i=0;i<volumeStructArraySize;i++) {
   
      if (volumeStructArray[i].volumePtr == volumePtr) { // Volume is already in the array!
      
         if ((volumeStructArray[i].n+i) == volumeStructArraySize) return R; // Sub-array is already at the end!

         memcpy(volumeStructArray+volumeStructArraySize,volumeStructArray+i,sizeof(volumeStruct)*volumeStructArray[i].n); // Copy sub-array to the end!
         volumeStructArraySize += volumeStructArray[i].n;
	 return R;
      }
   }

   xercesc::DOMElement* volumeElement = newElement("volume");
   volumeElement->setAttributeNode(newAttribute("name",volumePtr->GetName()));
   xercesc::DOMElement* materialrefElement = newElement("materialref");
   volumeElement->appendChild(materialrefElement);
   materialrefElement->setAttributeNode(newAttribute("ref",volumePtr->GetMaterial()->GetName()));
   xercesc::DOMElement* solidrefElement = newElement("solidref");
   volumeElement->appendChild(solidrefElement);
   solidrefElement->setAttributeNode(newAttribute("ref",solidPtr->GetName()));

   int sizeBefore = volumeStructArraySize;

   volumeStructArray[volumeStructArraySize].volumePtr = volumePtr;
   volumeStructArray[volumeStructArraySize].volumeElement = volumeElement;
   volumeStructArraySize++;

   G4Transform3D invR = R.inverse();

   const G4int daughterCount = volumePtr->GetNoDaughters();

   for (G4int i=0;i<daughterCount;i++) {
   
      const G4VPhysicalVolume* physvol = volumePtr->GetDaughter(i);
      physvolWrite(volumeElement,physvol,invR);
   }

   volumeStructArray[sizeBefore].n = volumeStructArraySize - sizeBefore;

   return R;
}

void G4GDMLWriteStructure::structureWrite(xercesc::DOMElement* gdmlElement,const G4LogicalVolume* worldvol) {

   G4cout << "Writing structure..." << G4endl;

   xercesc::DOMElement* structureElement = newElement("structure");
   gdmlElement->appendChild(structureElement);

   volumeStructArray = new volumeStruct[3000000];
   volumeStructArraySize = 0;

   if (volumeStructArray == 0) G4Exception("Not enough memory!");

   volumeWrite(worldvol);

   std::map<const G4LogicalVolume*,int> volumeMap;

   for (int i=0;i<volumeStructArraySize;i++) {
   
      int index = volumeStructArraySize-i-1;

      if (volumeMap.find(volumeStructArray[index].volumePtr) != volumeMap.end()) continue; // already printed!
   
      structureElement->appendChild(volumeStructArray[index].volumeElement);
   
      volumeMap[volumeStructArray[index].volumePtr] = 0;
   }

   delete [] volumeStructArray;
}
