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
// $Id: G4GDMLWriteStructure.cc,v 1.74.2.1 2009/08/11 08:27:49 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-02 $
//
// class G4GDMLWriteStructure Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLWriteStructure.hh"

G4GDMLWriteStructure::G4GDMLWriteStructure() : G4GDMLWriteParamvol()
{
}

G4GDMLWriteStructure::~G4GDMLWriteStructure()
{
}

void
G4GDMLWriteStructure::DivisionvolWrite(xercesc::DOMElement* volumeElement,
                                       const G4PVDivision* const divisionvol)
{
   EAxis axis = kUndefined;
   G4int number = 0;
   G4double width = 0.0;
   G4double offset = 0.0;
   G4bool consuming = false;

   divisionvol->GetReplicationData(axis,number,width,offset,consuming);

   G4String unitString("mm");
   G4String axisString("kUndefined");
   if (axis==kXAxis) { axisString = "kXAxis"; }
   else if (axis==kYAxis) { axisString = "kYAxis"; }
   else if (axis==kZAxis) { axisString = "kZAxis"; }
   else if (axis==kRho)   { axisString = "kRho";     }
   else if (axis==kPhi)   { axisString = "kPhi"; unitString = "degree"; }

   const G4String name
         = GenerateName(divisionvol->GetName(),divisionvol);
   const G4String volumeref
         = GenerateName(divisionvol->GetLogicalVolume()->GetName(),
                        divisionvol->GetLogicalVolume());

   xercesc::DOMElement* divisionvolElement = NewElement("divisionvol");
   divisionvolElement->setAttributeNode(NewAttribute("axis",axisString));
   divisionvolElement->setAttributeNode(NewAttribute("number",number));
   divisionvolElement->setAttributeNode(NewAttribute("width",width));
   divisionvolElement->setAttributeNode(NewAttribute("offset",offset));
   divisionvolElement->setAttributeNode(NewAttribute("unit",unitString));
   xercesc::DOMElement* volumerefElement = NewElement("volumeref");
   volumerefElement->setAttributeNode(NewAttribute("ref",volumeref));
   divisionvolElement->appendChild(volumerefElement);
   volumeElement->appendChild(divisionvolElement);
}

void G4GDMLWriteStructure::PhysvolWrite(xercesc::DOMElement* volumeElement,
                                        const G4VPhysicalVolume* const physvol,
                                        const G4Transform3D& T,
                                        const G4String& ModuleName)
{
   HepGeom::Scale3D scale;
   HepGeom::Rotate3D rotate;
   HepGeom::Translate3D translate;

   T.getDecomposition(scale,rotate,translate);

   const G4ThreeVector scl(scale(0,0),scale(1,1),scale(2,2));
   const G4ThreeVector rot = GetAngles(rotate.getRotation());
   const G4ThreeVector pos = T.getTranslation();

   const G4String name = GenerateName(physvol->GetName(),physvol);

   xercesc::DOMElement* physvolElement = NewElement("physvol");
   physvolElement->setAttributeNode(NewAttribute("name",name));
   volumeElement->appendChild(physvolElement);

   const G4String volumeref
         = GenerateName(physvol->GetLogicalVolume()->GetName(),
                        physvol->GetLogicalVolume());

   if (ModuleName.empty())
   {
      xercesc::DOMElement* volumerefElement = NewElement("volumeref");
      volumerefElement->setAttributeNode(NewAttribute("ref",volumeref));
      physvolElement->appendChild(volumerefElement);
   }
   else
   {
      xercesc::DOMElement* fileElement = NewElement("file");
      fileElement->setAttributeNode(NewAttribute("name",ModuleName));
      fileElement->setAttributeNode(NewAttribute("volname",volumeref));
      physvolElement->appendChild(fileElement);
   }

   if (std::fabs(pos.x()) > kLinearPrecision
    || std::fabs(pos.y()) > kLinearPrecision
    || std::fabs(pos.z()) > kLinearPrecision)
   {
     PositionWrite(physvolElement,name+"_pos",pos);
   }
   if (std::fabs(rot.x()) > kAngularPrecision
    || std::fabs(rot.y()) > kAngularPrecision
    || std::fabs(rot.z()) > kAngularPrecision)
   {
     RotationWrite(physvolElement,name+"_rot",rot);
   }
   if (std::fabs(scl.x()-1.0) > kRelativePrecision
    || std::fabs(scl.y()-1.0) > kRelativePrecision
    || std::fabs(scl.z()-1.0) > kRelativePrecision)
   {
     ScaleWrite(physvolElement,name+"_scl",scl);
   }
}

void G4GDMLWriteStructure::ReplicavolWrite(xercesc::DOMElement* volumeElement,
                                     const G4VPhysicalVolume* const replicavol)
{
   EAxis axis = kUndefined;
   G4int number = 0;
   G4double width = 0.0;
   G4double offset = 0.0;
   G4bool consuming = false;
   G4String unitString("mm");

   replicavol->GetReplicationData(axis,number,width,offset,consuming);

   const G4String volumeref
         = GenerateName(replicavol->GetLogicalVolume()->GetName(),
                        replicavol->GetLogicalVolume());

   xercesc::DOMElement* replicavolElement = NewElement("replicavol");
   replicavolElement->setAttributeNode(NewAttribute("number",number));
   xercesc::DOMElement* volumerefElement = NewElement("volumeref");
   volumerefElement->setAttributeNode(NewAttribute("ref",volumeref));
   replicavolElement->appendChild(volumerefElement);
   xercesc::DOMElement* replicateElement = NewElement("replicate_along_axis");
   replicavolElement->appendChild(replicateElement);

   xercesc::DOMElement* dirElement = NewElement("direction");
   if(axis==kXAxis)
     { dirElement->setAttributeNode(NewAttribute("x","1")); }
   else if(axis==kYAxis)
     { dirElement->setAttributeNode(NewAttribute("y","1")); }
   else if(axis==kZAxis)
     { dirElement->setAttributeNode(NewAttribute("z","1")); }
   else if(axis==kRho)
     { dirElement->setAttributeNode(NewAttribute("rho","1")); }
   else if(axis==kPhi)
     { dirElement->setAttributeNode(NewAttribute("phi","1")); }
   replicateElement->appendChild(dirElement);

   xercesc::DOMElement* widthElement = NewElement("width");
   widthElement->setAttributeNode(NewAttribute("value",width));
   widthElement->setAttributeNode(NewAttribute("unit",unitString));
   replicateElement->appendChild(widthElement);

   xercesc::DOMElement* offsetElement = NewElement("offset");
   offsetElement->setAttributeNode(NewAttribute("value",offset));
   offsetElement->setAttributeNode(NewAttribute("unit",unitString));
   replicateElement->appendChild(offsetElement);

   volumeElement->appendChild(replicavolElement);
}

void G4GDMLWriteStructure::StructureWrite(xercesc::DOMElement* gdmlElement)
{
   G4cout << "G4GDML: Writing structure..." << G4endl;

   structureElement = NewElement("structure");
   gdmlElement->appendChild(structureElement);
}

G4Transform3D G4GDMLWriteStructure::
TraverseVolumeTree(const G4LogicalVolume* const volumePtr, const G4int depth)
{
   if (VolumeMap().find(volumePtr) != VolumeMap().end())
   {
     return VolumeMap()[volumePtr]; // Volume is already processed
   }

   G4VSolid* solidPtr = volumePtr->GetSolid();
   G4Transform3D R,invR;
   G4int trans=0;

   while (true) // Solve possible displacement/reflection
   {            // of the referenced solid!
      if (trans>maxTransforms)
      {
        G4String ErrorMessage = "Referenced solid in volume '"
                              + volumePtr->GetName()
                              + "' was displaced/reflected too many times!";
        G4Exception("G4GDMLWriteStructure::TraverseVolumeTree()",
                    "InvalidSetup", FatalException, ErrorMessage);
      }

      if (G4ReflectedSolid* refl = dynamic_cast<G4ReflectedSolid*>(solidPtr))
      {
         R = R*refl->GetTransform3D();
         solidPtr = refl->GetConstituentMovedSolid();
         trans++;
         continue;
      }

      if (G4DisplacedSolid* disp = dynamic_cast<G4DisplacedSolid*>(solidPtr))
      {
         R = R*G4Transform3D(disp->GetObjectRotation(),
                             disp->GetObjectTranslation());
         solidPtr = disp->GetConstituentMovedSolid();
         trans++;
         continue;
      }

      break;
   }

   // Only compute the inverse when necessary!
   //
   if (trans>0) { invR = R.inverse(); }

   const G4String name
         = GenerateName(volumePtr->GetName(),volumePtr);
   const G4String materialref
         = GenerateName(volumePtr->GetMaterial()->GetName(),
                        volumePtr->GetMaterial());
   const G4String solidref
         = GenerateName(solidPtr->GetName(),solidPtr);

   xercesc::DOMElement* volumeElement = NewElement("volume");
   volumeElement->setAttributeNode(NewAttribute("name",name));
   xercesc::DOMElement* materialrefElement = NewElement("materialref");
   materialrefElement->setAttributeNode(NewAttribute("ref",materialref));
   volumeElement->appendChild(materialrefElement);
   xercesc::DOMElement* solidrefElement = NewElement("solidref");
   solidrefElement->setAttributeNode(NewAttribute("ref",solidref));
   volumeElement->appendChild(solidrefElement);

   const G4int daughterCount = volumePtr->GetNoDaughters();

   for (G4int i=0;i<daughterCount;i++)   // Traverse all the children!
   {
      const G4VPhysicalVolume* const physvol = volumePtr->GetDaughter(i);
      const G4String ModuleName = Modularize(physvol,depth);

      G4Transform3D daughterR;

      if (ModuleName.empty())   // Check if subtree requested to be 
      {                         // a separate module!
         daughterR = TraverseVolumeTree(physvol->GetLogicalVolume(),depth+1);
      }
      else
      {   
         G4GDMLWriteStructure writer;
         daughterR = writer.Write(ModuleName,physvol->GetLogicalVolume(),
                                  SchemaLocation,depth+1);
      }

      if (const G4PVDivision* const divisionvol
         = dynamic_cast<const G4PVDivision*>(physvol)) // Is it division?
      {
         if (!G4Transform3D::Identity.isNear(invR*daughterR,kRelativePrecision))
         {
            G4String ErrorMessage = "Division volume in '"
                                  + name
                                  + "' can not be related to reflected solid!";
            G4Exception("G4GDMLWriteStructure::TraverseVolumeTree()",
                        "InvalidSetup", FatalException, ErrorMessage);
         }
         DivisionvolWrite(volumeElement,divisionvol); 
      } else 
      if (physvol->IsParameterised())   // Is it a paramvol?
      {
         if (!G4Transform3D::Identity.isNear(invR*daughterR,kRelativePrecision))
         {
            G4String ErrorMessage = "Parameterised volume in '"
                                  + name
                                  + "' can not be related to reflected solid!";
            G4Exception("G4GDMLWriteStructure::TraverseVolumeTree()",
                        "InvalidSetup", FatalException, ErrorMessage);
         }
         ParamvolWrite(volumeElement,physvol);
      } else
      if (physvol->IsReplicated())   // Is it a replicavol?
      {
         if (!G4Transform3D::Identity.isNear(invR*daughterR,kRelativePrecision))
         {
            G4String ErrorMessage = "Replica volume in '"
                                  + name
                                  + "' can not be related to reflected solid!";
            G4Exception("G4GDMLWriteStructure::TraverseVolumeTree()",
                        "InvalidSetup", FatalException, ErrorMessage);
         }
         ReplicavolWrite(volumeElement,physvol); 
      }
      else   // Is it a physvol?
      {
         G4RotationMatrix rot;

         if (physvol->GetFrameRotation() != 0)
         {
           rot = *(physvol->GetFrameRotation());
         }
         G4Transform3D P(rot,physvol->GetObjectTranslation());
         PhysvolWrite(volumeElement,physvol,invR*P*daughterR,ModuleName);
      }
   }

   structureElement->appendChild(volumeElement);
     // Append the volume AFTER traversing the children so that
     // the order of volumes will be correct!

   VolumeMap()[volumePtr] = R;

   AddExtension(volumeElement, volumePtr);
     // Add any possible user defined extension attached to a volume

   AddMaterial(volumePtr->GetMaterial());
     // Add the involved materials and solids!

   AddSolid(solidPtr);

   return R;
}
