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

void G4GDMLWriteStructure::volumeWrite(const G4LogicalVolume* volumePtr) {

   if (volumeMap.find(volumePtr) != volumeMap.end()) { // Volume is already walked!
     
      volumeRecord volrec = volumeMap[volumePtr];
      volumeListType tempList(volrec.begin,volrec.end);
      volumeList.erase(volrec.begin,volrec.end);

      for (volumeListType::iterator i=tempList.begin();i != tempList.end();i++) {
      
         volumeList.push_back(*i);
      }

      return;
   }

   volumeRecord volrec;
   volrec.begin = volumeList.end();
   volrec.volumeElement = newElement("volume");
   volrec.volumeElement->setAttributeNode(newAttribute("name",volumePtr->GetName()));

   volumeList.push_back(volumePtr);

   const G4int daughterCount = volumePtr->GetNoDaughters();

   for (G4int i=0;i<daughterCount;i++) {
   
      const G4VPhysicalVolume* physvol = volumePtr->GetDaughter(i);
      const G4LogicalVolume* logvol = physvol->GetLogicalVolume();

      volumeWrite(logvol);
   }

   volrec.end = volumeList.end();
   volumeMap[volumePtr] = volrec;
}

void G4GDMLWriteStructure::structureWrite(xercesc::DOMElement* gdmlElement,const G4LogicalVolume* worldvol) {

   volumeList.reserve(G4LogicalVolumeStore::GetInstance()->size()); // IMPORTANT!!! There can be large datasets!

   xercesc::DOMElement* structureElement = newElement("structure");
   gdmlElement->appendChild(structureElement);

   volumeWrite(worldvol);

   for (volumeListType::reverse_iterator i=volumeList.rbegin();i != volumeList.rend();i++) {
   
      if (volumeMap.find(*i) == volumeMap.end()) G4Exception("GDML WRITER: Error in algorhythm! Volume in list is not in the map!");

      structureElement->appendChild(volumeMap[*i].volumeElement);
   }
}
