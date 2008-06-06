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
   divisionvolElement->setAttributeNode(newAttribute("name",divisionvol->GetName()));
   divisionvolElement->setAttributeNode(newAttribute("axis",axisString));
   divisionvolElement->setAttributeNode(newAttribute("number",number));
   divisionvolElement->setAttributeNode(newAttribute("width",width));
   divisionvolElement->setAttributeNode(newAttribute("offset",offset));
   divisionvolElement->setAttributeNode(newAttribute("unit",unitString));

   xercesc::DOMElement* volumerefElement = newElement("volumeref");
   divisionvolElement->appendChild(volumerefElement);
   volumerefElement->setAttributeNode(newAttribute("ref",divisionvol->GetLogicalVolume()->GetName()));
}

void G4GDMLWriteStructure::physvolWrite(xercesc::DOMElement* volumeElement,const G4VPhysicalVolume* const physvol,const G4Transform3D& T,const G4String& ModuleName) {

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

   if (ModuleName.empty()) {

      xercesc::DOMElement* volumerefElement = newElement("volumeref");
      volumerefElement->setAttributeNode(newAttribute("ref",physvol->GetLogicalVolume()->GetName()));
      physvolElement->appendChild(volumerefElement);
   } else {

      xercesc::DOMElement* fileElement = newElement("file");
      fileElement->setAttributeNode(newAttribute("name",ModuleName));
      fileElement->setAttributeNode(newAttribute("volname",physvol->GetLogicalVolume()->GetName()));
      physvolElement->appendChild(fileElement);
   }

   if (fabs(scl.x()-1.0) > kRelativePrecision || fabs(scl.y()-1.0) > kRelativePrecision || fabs(scl.z()-1.0) > kRelativePrecision) scaleWrite(physvolElement,"",scl);
   if (fabs(rot.x()) > kAngularPrecision || fabs(rot.y()) > kAngularPrecision || fabs(rot.z()) > kAngularPrecision) rotationWrite(physvolElement,"",rot);
   if (fabs(pos.x()) > kLinearPrecision || fabs(pos.y()) > kLinearPrecision || fabs(pos.z()) > kLinearPrecision) positionWrite(physvolElement,"",pos);
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
   replicavolElement->setAttributeNode(newAttribute("name",replicavol->GetName()));
   replicavolElement->setAttributeNode(newAttribute("axis",axisString));
   replicavolElement->setAttributeNode(newAttribute("number",number));
   replicavolElement->setAttributeNode(newAttribute("width",width));
   replicavolElement->setAttributeNode(newAttribute("offset",offset));
   replicavolElement->setAttributeNode(newAttribute("unit",unitString));

   xercesc::DOMElement* volumerefElement = newElement("volumeref");
   replicavolElement->appendChild(volumerefElement);
   volumerefElement->setAttributeNode(newAttribute("ref",replicavol->GetLogicalVolume()->GetName()));
}

void G4GDMLWriteStructure::structureWrite(xercesc::DOMElement* gdmlElement) {

   G4cout << "Writing structure..." << G4endl;

   structureElement = newElement("structure");
   gdmlElement->appendChild(structureElement);
}

G4Transform3D G4GDMLWriteStructure::TraverseVolumeTree(const G4LogicalVolume* const volumePtr) {

   if (volumeMap.find(volumePtr) != volumeMap.end()) return volumeMap[volumePtr]; // Volume is already processed

   G4Transform3D R;
   int displaced = 0;
   G4VSolid* solidPtr = volumePtr->GetSolid();

   while (true) { // Solve possible displacement/reflection of the referenced solid!
   
      if (displaced>maxDisplacements) G4Exception("GDML Writer: ERROR! Referenced solid in volume '"+volumePtr->GetName()+"' was displaced/reflected too many times!");
   
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

   xercesc::DOMElement* volumeElement = newElement("volume");
   volumeElement->setAttributeNode(newAttribute("name",volumePtr->GetName()));
   xercesc::DOMElement* materialrefElement = newElement("materialref");
   volumeElement->appendChild(materialrefElement);
   materialrefElement->setAttributeNode(newAttribute("ref",volumePtr->GetMaterial()->GetName()));
   xercesc::DOMElement* solidrefElement = newElement("solidref");
   volumeElement->appendChild(solidrefElement);
   solidrefElement->setAttributeNode(newAttribute("ref",solidPtr->GetName()));

   G4Transform3D invR = R.inverse();

   const G4int daughterCount = volumePtr->GetNoDaughters();

   for (G4int i=0;i<daughterCount;i++) { // Traverse all the children!
   
      const G4VPhysicalVolume* const physvol = volumePtr->GetDaughter(i);
      G4String ModuleName = G4GDMLWrite::GetModule(physvol); 
      G4Transform3D daughterR;
      
      if (ModuleName.empty()) { // Check if the subtree starting with this volume is requested to be a spearate module or not!
      
         daughterR = TraverseVolumeTree(physvol->GetLogicalVolume());
      } else {

         G4GDMLWriteStructure writer;
         daughterR = writer.Write(ModuleName,physvol->GetLogicalVolume());
      }
      
      if (const G4PVDivision* const divisionvol = dynamic_cast<const G4PVDivision* const>(physvol)) { 
      
         if (!G4Transform3D::Identity.isNear(invR*daughterR,kRelativePrecision)) 
	    G4Exception("GDML Writer: Error! divisionvol in '"+volumePtr->GetName()+"' can not be related to reflected solid!");

         divisionvolWrite(volumeElement,divisionvol); 
      } else 
      if (physvol->IsParameterised()) { 
       
         if (!G4Transform3D::Identity.isNear(invR*daughterR,kRelativePrecision)) 
	    G4Exception("GDML Writer: Error! paramvol in '"+volumePtr->GetName()+"' can not be related to reflected solid!");

         paramvolWrite(volumeElement,physvol);
      } else
      if (physvol->IsReplicated()) { 

         if (!G4Transform3D::Identity.isNear(invR*daughterR,kRelativePrecision))
	    G4Exception("GDML Writer: Error! replicavol in '"+volumePtr->GetName()+"' can not be related to reflected solid!");

         replicavolWrite(volumeElement,physvol); 
      } else {
   
         G4Transform3D P(physvol->GetObjectRotationValue().inverse(),physvol->GetObjectTranslation());
         physvolWrite(volumeElement,physvol,invR*P*daughterR,ModuleName);
      }
   }

   structureElement->appendChild(volumeElement); // Append the volume AFTER traversing the children so that the order of volumes will be correct!
   volumeMap[volumePtr] = R;

   G4GDMLWriteMaterials::AddMaterial(volumePtr->GetMaterial());   // Only the involved materials and solids are added!
   G4GDMLWriteSolids::AddSolid(solidPtr);

   return R;
}
