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

void G4GDMLWriteStructure::physvolWrite(xercesc::DOMElement* volumeElement,const G4VPhysicalVolume* const physvol) {

   G4ThreeVector scl(1.0,1.0,1.0);
   G4ThreeVector rot = getAngles(physvol->GetObjectRotationValue().inverse());
   G4ThreeVector pos = physvol->GetObjectTranslation();

   G4LogicalVolume* logvol = physvol->GetLogicalVolume();
   G4VSolid* solidPtr = logvol->GetSolid();
   G4String volumeref = logvol->GetName();

   volumeWrite(logvol);   

   if (const G4ReflectedSolid* refl = dynamic_cast<const G4ReflectedSolid*>(solidPtr)) {

      volumeref.remove(volumeref.length()-5,5); // Remove "_refl"
      G4Transform3D scale = refl->GetTransform3D();
      	 
      scl.setX(scale(0,0)); // We assume that this is a pure scaling transformation!!!
      scl.setY(scale(1,1));
      scl.setZ(scale(2,2));
   }

   xercesc::DOMElement* physvolElement = newElement("physvol");
   volumeElement->appendChild(physvolElement);
   physvolElement->setAttributeNode(newAttribute("name",physvol->GetName()));

   xercesc::DOMElement* volumerefElement = newElement("volumeref");
   physvolElement->appendChild(volumerefElement);
   volumerefElement->setAttributeNode(newAttribute("ref",volumeref));

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

void G4GDMLWriteStructure::volumeWrite(const G4LogicalVolume* const volumePtr) {

   G4VSolid* solidPtr = volumePtr->GetSolid();

   if (dynamic_cast<const G4ReflectedSolid*>(solidPtr)) return; // Reflected solid is replaced with scale transformation!
      
   xercesc::DOMElement* volumeElement = newElement("volume");

   volumeElementPair pair;
   pair.key = const_cast<G4LogicalVolume*>(volumePtr);
   pair.value = volumeElement;
   volumeElementList.push_back(pair);

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
/*   
      if (const G4PVDivision* const divisionvol = dynamic_cast<const G4PVDivision*>(physvol)) divisionvolWrite(volumeElement,divisionvol); else
      if (physvol->IsParameterised()) paramvolWrite(volumeElement,physvol); else
      if (physvol->IsReplicated()) replicavolWrite(volumeElement,physvol); else*/
      physvolWrite(volumeElement,physvol);
   }
}

void G4GDMLWriteStructure::structureWrite(xercesc::DOMElement* gdmlElement,const G4LogicalVolume* const worldvol) {

   xercesc::DOMElement* structureElement = newElement("structure");
   gdmlElement->appendChild(structureElement);

   volumeWrite(worldvol);

   for (int i=volumeElementList.size()-1;i>=0;i--) // Write back the volumes in reverse order!
      structureElement->appendChild(volumeElementList[i].value);  
}
