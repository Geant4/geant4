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

   const G4String name = GenerateName(divisionvol->GetName(),divisionvol);
   const G4String volumeref = GenerateName(divisionvol->GetLogicalVolume()->GetName(),divisionvol->GetLogicalVolume());

   xercesc::DOMElement* divisionvolElement = newElement("divisionvol");
   divisionvolElement->setAttributeNode(newAttribute("name",name));
   divisionvolElement->setAttributeNode(newAttribute("axis",axisString));
   divisionvolElement->setAttributeNode(newAttribute("number",number));
   divisionvolElement->setAttributeNode(newAttribute("width",width));
   divisionvolElement->setAttributeNode(newAttribute("offset",offset));
   divisionvolElement->setAttributeNode(newAttribute("unit",unitString));
   xercesc::DOMElement* volumerefElement = newElement("volumeref");
   volumerefElement->setAttributeNode(newAttribute("ref",volumeref));
   divisionvolElement->appendChild(volumerefElement);
   volumeElement->appendChild(divisionvolElement);
}

void G4GDMLWriteStructure::physvolWrite(xercesc::DOMElement* volumeElement,const G4VPhysicalVolume* const physvol,const G4Transform3D& T,const G4String& ModuleName) {

   HepGeom::Scale3D scale;
   HepGeom::Rotate3D rotate;
   HepGeom::Translate3D translate;

   T.getDecomposition(scale,rotate,translate);

   const G4ThreeVector scl(scale(0,0),scale(1,1),scale(2,2));
   const G4ThreeVector rot = getAngles(rotate.getRotation());
   const G4ThreeVector pos = T.getTranslation();

   const G4String name = GenerateName(physvol->GetName(),physvol);

   xercesc::DOMElement* physvolElement = newElement("physvol");
   physvolElement->setAttributeNode(newAttribute("name",name));
   volumeElement->appendChild(physvolElement);

   const G4String volumeref = GenerateName(physvol->GetLogicalVolume()->GetName(),physvol->GetLogicalVolume());

   if (ModuleName.empty()) {

      xercesc::DOMElement* volumerefElement = newElement("volumeref");
      volumerefElement->setAttributeNode(newAttribute("ref",volumeref));
      physvolElement->appendChild(volumerefElement);
   } else {

      xercesc::DOMElement* fileElement = newElement("file");
      fileElement->setAttributeNode(newAttribute("name",ModuleName));
      fileElement->setAttributeNode(newAttribute("volname",volumeref));
      physvolElement->appendChild(fileElement);
   }

   std::stringstream ptr_stream;
   ptr_stream << physvol;
   const G4String ptr_string = ptr_stream.str();

   if (fabs(pos.x()) > kLinearPrecision || fabs(pos.y()) > kLinearPrecision || fabs(pos.z()) > kLinearPrecision) positionWrite(physvolElement,"position"+ptr_string,pos);
   if (fabs(rot.x()) > kAngularPrecision || fabs(rot.y()) > kAngularPrecision || fabs(rot.z()) > kAngularPrecision) rotationWrite(physvolElement,"rotation"+ptr_string,rot);
   if (fabs(scl.x()-1.0) > kRelativePrecision || fabs(scl.y()-1.0) > kRelativePrecision || fabs(scl.z()-1.0) > kRelativePrecision) scaleWrite(physvolElement,"scale"+ptr_string,scl);
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

   const G4String name = GenerateName(replicavol->GetName(),replicavol);
   const G4String volumeref = GenerateName(replicavol->GetLogicalVolume()->GetName(),replicavol->GetLogicalVolume());

   xercesc::DOMElement* replicavolElement = newElement("replicavol");
   replicavolElement->setAttributeNode(newAttribute("name",name));
   replicavolElement->setAttributeNode(newAttribute("axis",axisString));
   replicavolElement->setAttributeNode(newAttribute("number",number));
   replicavolElement->setAttributeNode(newAttribute("width",width));
   replicavolElement->setAttributeNode(newAttribute("offset",offset));
   replicavolElement->setAttributeNode(newAttribute("unit",unitString));
   xercesc::DOMElement* volumerefElement = newElement("volumeref");
   volumerefElement->setAttributeNode(newAttribute("ref",volumeref));
   replicavolElement->appendChild(volumerefElement);
   volumeElement->appendChild(replicavolElement);
}

void G4GDMLWriteStructure::structureWrite(xercesc::DOMElement* gdmlElement) {

   G4cout << "G4GDML: Writing structure..." << G4endl;

   structureElement = newElement("structure");
   gdmlElement->appendChild(structureElement);
}

G4Transform3D G4GDMLWriteStructure::TraverseVolumeTree(const G4LogicalVolume* const volumePtr,const G4int depth) {

   if (volumeMap.find(volumePtr) != volumeMap.end()) return volumeMap[volumePtr]; // Volume is already processed

   G4VSolid* solidPtr = volumePtr->GetSolid();
   G4Transform3D R,invR;
   G4int reflected = 0;

   while (true) { // Solve possible displacement/reflection of the referenced solid!
   
      if (reflected>maxReflections) G4Exception("G4GDML: ERROR! Referenced solid in volume '"+volumePtr->GetName()+"' was displaced/reflected too many times!");
   
      if (G4ReflectedSolid* refl = dynamic_cast<G4ReflectedSolid*>(solidPtr)) {
   
         R = R*refl->GetTransform3D();
         solidPtr = refl->GetConstituentMovedSolid();
         reflected++;
         continue;
      }

      if (G4DisplacedSolid* disp = dynamic_cast<G4DisplacedSolid*>(solidPtr)) {
      
         R = R*G4Transform3D(disp->GetObjectRotation(),disp->GetObjectTranslation());
         solidPtr = disp->GetConstituentMovedSolid();
         reflected++;
         continue;
      }

      break;
   }

   if (reflected>0) invR = R.inverse(); // Only compute the inverse when necessary!

   const G4String name = GenerateName(volumePtr->GetName(),volumePtr);
   const G4String materialref = GenerateName(volumePtr->GetMaterial()->GetName(),volumePtr->GetMaterial());
   const G4String solidref = GenerateName(solidPtr->GetName(),solidPtr);

   xercesc::DOMElement* volumeElement = newElement("volume");
   volumeElement->setAttributeNode(newAttribute("name",name));
   xercesc::DOMElement* materialrefElement = newElement("materialref");
   materialrefElement->setAttributeNode(newAttribute("ref",materialref));
   volumeElement->appendChild(materialrefElement);
   xercesc::DOMElement* solidrefElement = newElement("solidref");
   solidrefElement->setAttributeNode(newAttribute("ref",solidref));
   volumeElement->appendChild(solidrefElement);

   const G4int daughterCount = volumePtr->GetNoDaughters();

   for (G4int i=0;i<daughterCount;i++) { // Traverse all the children!
   
      const G4VPhysicalVolume* const physvol = volumePtr->GetDaughter(i);
      const G4String ModuleName = Modularize(physvol,depth);
      G4Transform3D daughterR;
      
      if (ModuleName.empty()) { // Check if subtree requested to be a separate module!

         daughterR = TraverseVolumeTree(physvol->GetLogicalVolume(),depth+1);
      } else {
         
	 G4GDMLWriteStructure writer;
         daughterR = writer.Write(ModuleName,physvol->GetLogicalVolume(),depth+1);
      }

      if (const G4PVDivision* const divisionvol = dynamic_cast<const G4PVDivision* const>(physvol)) { // Is it a divisionvol?
      
         if (!G4Transform3D::Identity.isNear(invR*daughterR,kRelativePrecision)) 
	    G4Exception("G4GDML: ERROR! divisionvol in '"+name+"' can not be related to reflected solid!");

         divisionvolWrite(volumeElement,divisionvol); 
      } else 
      if (physvol->IsParameterised()) { // Is it a paramvol?
       
         if (!G4Transform3D::Identity.isNear(invR*daughterR,kRelativePrecision)) 
	    G4Exception("G4GDML: ERROR! paramvol in '"+name+"' can not be related to reflected solid!");

         paramvolWrite(volumeElement,physvol);
      } else
      if (physvol->IsReplicated()) { // Is it a replicavol?

         if (!G4Transform3D::Identity.isNear(invR*daughterR,kRelativePrecision))
	    G4Exception("G4GDML: ERROR! replicavol in '"+name+"' can not be related to reflected solid!");

         replicavolWrite(volumeElement,physvol); 
      } else { // Is it a physvol?
   
         G4RotationMatrix rot;

         if (physvol->GetFrameRotation() != 0) rot = *(physvol->GetFrameRotation());
   
         G4Transform3D P(rot,physvol->GetObjectTranslation());
         physvolWrite(volumeElement,physvol,invR*P*daughterR,ModuleName);
      }
   }

   structureElement->appendChild(volumeElement); // Append the volume AFTER traversing the children so that the order of volumes will be correct!

   volumeMap[volumePtr] = R;

   G4GDMLWriteMaterials::AddMaterial(volumePtr->GetMaterial());   // Add the involved materials and solids!
   G4GDMLWriteSolids::AddSolid(solidPtr);

   return R;
}
