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

void G4GDMLWriteStructure::replicavolWrite(xercesc::DOMElement* volumeElement,const G4VPhysicalVolume* const replicavol) {

   EAxis axis = kUndefined;
   G4int number = 0;
   G4double width = 0.0;
   G4double offset = 0.0;
   G4bool consuming = false;

   replicavol->GetReplicationData(axis,number,width,offset,consuming);

   G4String unitString("mm");
   G4String axisString("kUndefined");
   if (axis==kXAxis) axisString = "kXAxis"; else
   if (axis==kYAxis) axisString = "kYAxis"; else
   if (axis==kZAxis) axisString = "kZAxis"; else
   if (axis==kRho) { axisString = "kRho"; unitString = "degree"; } else
   if (axis==kPhi) { axisString = "kPhi"; unitString = "degree"; }

   xercesc::DOMElement* replicavolElement = newElement("replicavol");
   volumeElement->appendChild(replicavolElement);
   replicavolElement->setAttributeNode(newAttribute("axis",axisString));
   replicavolElement->setAttributeNode(newAttribute("number",number));
   replicavolElement->setAttributeNode(newAttribute("width",width));
   replicavolElement->setAttributeNode(newAttribute("offset",offset));
   replicavolElement->setAttributeNode(newAttribute("unit",unitString));

   xercesc::DOMElement* volumerefElement = newElement("volumeref");
   replicavolElement->appendChild(volumerefElement);
   volumerefElement->setAttributeNode(newAttribute("ref",replicavol->GetLogicalVolume()->GetName()));
}

void G4GDMLWriteStructure::divisionvolWrite(xercesc::DOMElement* volumeElement,const G4PVDivision* const divisionvol) {

   EAxis axis = kUndefined;
   G4int number = 0;
   G4double width = 0.0;
   G4double offset = 0.0;
   G4bool consuming = false;

   divisionvol->GetReplicationData(axis,number,width,offset,consuming);

   G4String unitString("mm");
   G4String axisString("kUndefined");
   if (axis==kXAxis) axisString = "kXAxis"; else
   if (axis==kYAxis) axisString = "kYAxis"; else
   if (axis==kZAxis) axisString = "kZAxis"; else
   if (axis==kRho) { axisString = "kRho"; unitString = "degree"; } else
   if (axis==kPhi) { axisString = "kPhi"; unitString = "degree"; }

   xercesc::DOMElement* divisionvolElement = newElement("divisionvol");
   volumeElement->appendChild(divisionvolElement);
   divisionvolElement->setAttributeNode(newAttribute("axis",axisString));
   divisionvolElement->setAttributeNode(newAttribute("number",number));
   divisionvolElement->setAttributeNode(newAttribute("width",width));
   divisionvolElement->setAttributeNode(newAttribute("offset",offset));
   divisionvolElement->setAttributeNode(newAttribute("unit",unitString));

   xercesc::DOMElement* volumerefElement = newElement("volumeref");
   divisionvolElement->appendChild(volumerefElement);
   volumerefElement->setAttributeNode(newAttribute("ref",divisionvol->GetLogicalVolume()->GetName()));
}

G4Transform3D G4GDMLWriteStructure::volumeWrite(const G4LogicalVolume* const volumePtr) {

   G4cout << "Walking volume: " << volumePtr->GetName() << G4endl;

   G4Transform3D R;

   const G4VSolid* solidPtr = volumePtr->GetSolid();
   if (const G4ReflectedSolid* refl = dynamic_cast<const G4ReflectedSolid*>(solidPtr)) {
   
      solidPtr = refl->GetConstituentMovedSolid();
      R = refl->GetTransform3D();
   }

   for (volumePtrListType::iterator i=volumePtrList.begin();i != volumePtrList.end();i++) {
   
      if (*i == volumePtr) {
      
         volumePtrList.erase(i);
         break;
      }
   }

   xercesc::DOMElement* volumeElement = newElement("volume");

//   if (volumePtrMap.find(volumePtr) != volumePtrMap.end()) G4Exception("ERROR! Volume is already added to map!"); 

   volumePtrList.insert(volumePtrList.begin(),volumePtr);
   volumePtrMap[volumePtr] = volumeElement;

   volumeElement->setAttributeNode(newAttribute("name",volumePtr->GetName()));

   xercesc::DOMElement* materialrefElement = newElement("materialref");
   materialrefElement->setAttributeNode(newAttribute("ref",volumePtr->GetMaterial()->GetName()));
   volumeElement->appendChild(materialrefElement);

   xercesc::DOMElement* solidrefElement = newElement("solidref");
   solidrefElement->setAttributeNode(newAttribute("ref",solidPtr->GetName()));
   volumeElement->appendChild(solidrefElement);

   const G4int daughterCount = volumePtr->GetNoDaughters();

   for (G4int i=0;i<daughterCount;i++) {
   
      const G4VPhysicalVolume* const physvol = volumePtr->GetDaughter(i);
      physvolWrite(volumeElement,physvol,R.inverse());
   }

   return R;
}

void G4GDMLWriteStructure::structureWrite(xercesc::DOMElement* gdmlElement,const G4LogicalVolume* const worldvol) {

   xercesc::DOMElement* structureElement = newElement("structure");
   gdmlElement->appendChild(structureElement);

   volumeWrite(worldvol);
 
   for (volumePtrListType::iterator i=volumePtrList.begin();i != volumePtrList.end();i++) {

      if (volumePtrMap.find(*i) == volumePtrMap.end()) G4Exception("ERROR! Volume is not added to map!"); 
      structureElement->appendChild(volumePtrMap[*i]);
   }
}
