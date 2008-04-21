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

void G4GDMLWriteStructure::physvolWrite(xercesc::DOMElement* volumeElement,const G4VPhysicalVolume* const physvol,const G4Transform3D& invR,const G4Transform3D& R) {

   G4Transform3D P(physvol->GetObjectRotationValue().inverse(),physvol->GetObjectTranslation());
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

G4Transform3D G4GDMLWriteStructure::volumeWrite(const G4LogicalVolume* volumePtr) {

   G4Transform3D R;
   G4VSolid* solidPtr = volumePtr->GetSolid();

   int displaced = 0;

   while (true) { // Solve possible displacement/reflection of the referenced solid!
   
      if (displaced>4) G4Exception("ERROR! Referenced solid in volume '"+volumePtr->GetName()+"' was displaced/reflected too many times!");
   
      if (G4ReflectedSolid* refl = dynamic_cast<G4ReflectedSolid*>(solidPtr)) {
   
         R = R*refl->GetTransform3D();
         solidPtr = refl->GetConstituentMovedSolid();
         displaced++;
         continue;
      }

      if (G4DisplacedSolid* disp = dynamic_cast<G4DisplacedSolid*>(solidPtr)) {
      
         R = R*G4Transform3D(disp->GetObjectRotation(),disp->GetObjectTranslation());
         solidPtr = disp->GetConstituentMovedSolid();
         displaced++;
         continue;
      }

      break;
   }
   
   for (int i=0;i<volumeArraySize;i++) {
   
      if (volumeArray[i]->volumePtr == volumePtr) { // Volume is already in the array!
      
         if ((volumeArray[i]->n+i) == volumeArraySize) return R; // Sub-array is already at the end!

         if ((volumeArraySize+volumeArray[i]->n) >= volumeArrayMaxSize) 
	    G4Exception("Error at sorting volumes! Volume array size is too small!");

         memcpy(volumeArray+volumeArraySize,volumeArray+i,sizeof(volumeStruct*)*volumeArray[i]->n); // Copy sub-array to the end!
         volumeArraySize += volumeArray[i]->n;
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

   volumeStruct* vols = new volumeStruct;
   vols->volumePtr = volumePtr;
   vols->volumeElement = volumeElement;
   vols->n = volumeArraySize;

   if (volumeArraySize >= volumeArrayMaxSize) G4Exception("Error at sorting volumes! Volume array size is too small!");

   volumeArray[volumeArraySize++] = vols;

   G4Transform3D invR = R.inverse();

   const G4int daughterCount = volumePtr->GetNoDaughters();

   for (G4int i=0;i<daughterCount;i++) {
   
      const G4VPhysicalVolume* physvol = volumePtr->GetDaughter(i);
      const G4Transform3D daughterR = volumeWrite(physvol->GetLogicalVolume());

      if (const G4PVDivision* const divisionvol = dynamic_cast<const G4PVDivision* const>(physvol)) { 
      
         if (!G4Transform3D::Identity.isNear(invR*daughterR)) G4Exception("Error! divisionvol in '"+volumePtr->GetName()+"' can not be related to reflected solid!");
         divisionvolWrite(volumeElement,divisionvol); 
      } else 
      if (physvol->IsParameterised()) { 
       
         if (!G4Transform3D::Identity.isNear(invR*daughterR)) G4Exception("Error! paramvol in '"+volumePtr->GetName()+"' can not be related to reflected solid!");
         paramvolWrite(volumeElement,physvol);
      } else
      if (physvol->IsReplicated()) { 

         if (!G4Transform3D::Identity.isNear(invR*daughterR)) G4Exception("Error! replicavol in '"+volumePtr->GetName()+"' can not be related to reflected solid!");
         replicavolWrite(volumeElement,physvol); 
      } else
      physvolWrite(volumeElement,physvol,invR,daughterR);
   }

   vols->n = volumeArraySize - vols->n; // Growth of array size after walking a subtree of a node = number of descendants of that node

   return R;
}

void G4GDMLWriteStructure::structureWrite(xercesc::DOMElement* gdmlElement,const G4LogicalVolume* worldvol) {

   G4cout << "Writing structure..." << G4endl;

   xercesc::DOMElement* structureElement = newElement("structure");
   gdmlElement->appendChild(structureElement);

   volumeArray = new volumeStruct*[volumeArrayMaxSize];
   volumeArraySize = 0;

   if (volumeArray == 0) G4Exception("Not enough memory for sorting volumes!");

   volumeWrite(worldvol);

   for (int i=volumeArraySize-1;i>=0;i--) { // Append all XML elements, representing volumes, in the correct (reverse) order!
   
      if (volumeArray[i]->volumePtr == 0) continue; // Append only once!

      structureElement->appendChild(volumeArray[i]->volumeElement);
      volumeArray[i]->volumePtr = 0;
   }

   delete [] volumeArray;
}
