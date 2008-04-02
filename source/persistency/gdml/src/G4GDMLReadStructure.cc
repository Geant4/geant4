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
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLReadStructure.hh"

void G4GDMLReadStructure::GeneratePhysvolName(G4VPhysicalVolume* physvol) {

   std::stringstream stream;
   stream << physvol->GetLogicalVolume()->GetName() << "_" << physvol;
   physvol->SetName(GenerateName(stream.str()));
}

void G4GDMLReadStructure::assemblyRead(const xercesc::DOMElement* const assemblyElement) {

   G4String name;

   const xercesc::DOMNamedNodeMap* const attributes = assemblyElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="name") name = GenerateName(attValue);
   }

   pAssembly = new G4AssemblyVolume();

   assemblyMap[name] = pAssembly;
}

G4GDMLAuxPairType G4GDMLReadStructure::auxiliaryRead(const xercesc::DOMElement* const auxiliaryElement) {

   G4GDMLAuxPairType auxpair;

   const xercesc::DOMNamedNodeMap* const attributes = auxiliaryElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="auxtype") auxpair.type = attValue; else
      if (attName=="auxvalue") auxpair.value = eval.Evaluate(attValue);
   }

   return auxpair;
}

void G4GDMLReadStructure::bordersurfaceRead(const xercesc::DOMElement* const bordersurfaceElement) {

   G4String name;
   G4VPhysicalVolume* pv1 = 0;
   G4VPhysicalVolume* pv2 = 0;
   G4SurfaceProperty* prop = 0;
   G4int index = 0;

   const xercesc::DOMNamedNodeMap* const attributes = bordersurfaceElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="name") name = GenerateName(attValue); else
      if (attName=="surfaceproperty") prop = getSurfaceProperty(GenerateName(attValue));
   }

   for (xercesc::DOMNode* iter = bordersurfaceElement->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag != "physvolref") continue; 
      
      if (index==0) { pv1 = getPhysvol(GenerateName(refRead(child))); index++; } else
      if (index==1) { pv2 = getPhysvol(GenerateName(refRead(child))); index++; } else
      break;
   }

   new G4LogicalBorderSurface(name,pv1,pv2,prop);
}

void G4GDMLReadStructure::divisionvolRead(const xercesc::DOMElement* const divisionvolElement) {

   G4double unit = 1.0;
   G4double width = 0.0;
   G4double offset = 0.0;
   G4int number = 0;
   EAxis axis = kUndefined;
   G4LogicalVolume* logvol = 0;
   
   const xercesc::DOMNamedNodeMap* const attributes = divisionvolElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="unit") unit = eval.Evaluate(attValue); else
      if (attName=="width") width = eval.Evaluate(attValue); else
      if (attName=="offset") offset = eval.Evaluate(attValue); else
      if (attName=="number") number = eval.EvaluateInteger(attValue); else
      if (attName=="axis") {
      
         if (attValue=="kXAxis") axis = kXAxis; else
         if (attValue=="kYAxis") axis = kYAxis; else
         if (attValue=="kZAxis") axis = kZAxis; else
         if (attValue=="kRho") axis = kRho; else
         if (attValue=="kPhi") axis = kPhi;
      }
   }

   width *= unit;
   offset *= unit;

   for (xercesc::DOMNode* iter = divisionvolElement->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="volumeref") logvol = getVolume(GenerateName(refRead(child)));
   }

   G4PVDivision* divisionvolPtr = 0;
   if (number != 0 && width == 0.0) divisionvolPtr = new G4PVDivision("",logvol,pMotherLogical,axis,number,offset); else
   if (number == 0 && width != 0.0) divisionvolPtr = new G4PVDivision("",logvol,pMotherLogical,axis,width,offset); else
   divisionvolPtr = new G4PVDivision("",logvol,pMotherLogical,axis,number,width,offset); 

   GeneratePhysvolName(divisionvolPtr);
}

G4LogicalVolume* G4GDMLReadStructure::fileRead(const xercesc::DOMElement* const fileElement) {

   G4String name;
   G4String volname;

   const xercesc::DOMNamedNodeMap* const attributes = fileElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="name") name = attValue; else
      if (attName=="volname") volname = attValue;
   }

   G4GDMLReadStructure structure; // We create a new structure with a new evaluator
   
   structure.Read(name,true); // true: it is an external file

   if (volname.empty()) return structure.getVolume(structure.getSetup("Default"));
   else return structure.getVolume(structure.GenerateName(volname));
}

void G4GDMLReadStructure::physvolRead(const xercesc::DOMElement* const physvolElement) {

   G4String name;
   G4LogicalVolume* logvol = 0;
   G4ThreeVector position(0.0,0.0,0.0);
   G4ThreeVector rotation(0.0,0.0,0.0);
   G4ThreeVector scale(1.0,1.0,1.0);

   const xercesc::DOMNamedNodeMap* const attributes = physvolElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="name") name = attValue;
   }

   for (xercesc::DOMNode* iter = physvolElement->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="file") logvol = fileRead(child); else
      if (tag=="volumeref") logvol = getVolume(GenerateName(refRead(child))); else
      if (tag=="position") vectorRead(child,position); else
      if (tag=="rotation") vectorRead(child,rotation); else
      if (tag=="scale") vectorRead(child,scale); else
      if (tag=="positionref") position = getPosition(GenerateName(refRead(child))); else
      if (tag=="rotationref") rotation = getRotation(GenerateName(refRead(child))); else
      if (tag=="scaleref") scale = getScale(GenerateName(refRead(child))); else
      G4Exception("GDML Reader: ERROR! Unknown tag in physvol: "+tag);
   }

   G4Transform3D transform(getRotationMatrix(rotation).inverse(),position);
   transform = transform*G4Scale3D(scale.x(),scale.y(),scale.z());

   G4PhysicalVolumesPair pair = G4ReflectionFactory::Instance()->Place(transform,GenerateName(name),logvol,pMotherLogical,false,0,false);

   if (name.empty()) {

      if (pair.first != 0) GeneratePhysvolName(pair.first);
      if (pair.second != 0) GeneratePhysvolName(pair.second);
   }
}

void G4GDMLReadStructure::replicavolRead(const xercesc::DOMElement* const replicavolElement) {

   G4double unit = 1.0;
   G4double width = 0.0;
   G4double offset = 0.0;
   G4int number = 0;
   EAxis axis = kUndefined;
   G4LogicalVolume* logvol = 0;

   const xercesc::DOMNamedNodeMap* const attributes = replicavolElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="unit") unit = eval.Evaluate(attValue); else
      if (attName=="width") width = eval.Evaluate(attValue); else
      if (attName=="offset") offset = eval.Evaluate(attValue); else
      if (attName=="number") number = eval.EvaluateInteger(attValue); else
      if (attName=="axis") {
      
         if (attValue=="kXAxis") axis = kXAxis; else
         if (attValue=="kYAxis") axis = kYAxis; else
         if (attValue=="kZAxis") axis = kZAxis; else
         if (attValue=="kRho") axis = kRho; else
         if (attValue=="kPhi") axis = kPhi;
      }
   }

   width *= unit;
   offset *= unit;

   for (xercesc::DOMNode* iter = replicavolElement->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="volumeref") logvol = getVolume(GenerateName(refRead(child)));
   }

   G4PVReplica* replicaPtr = new G4PVReplica("",logvol,pMotherLogical,axis,number,width,offset);
   GeneratePhysvolName(replicaPtr);
}

void G4GDMLReadStructure::volumeRead(const xercesc::DOMElement* const volumeElement) {

   G4VSolid* solidPtr = 0;
   G4Material* materialPtr = 0;
   G4GDMLAuxListType auxList;
   
   XMLCh *name_attr = xercesc::XMLString::transcode("name");
   G4String name = xercesc::XMLString::transcode(volumeElement->getAttribute(name_attr));
   xercesc::XMLString::release(&name_attr);

   for (xercesc::DOMNode* iter = volumeElement->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="auxiliary") auxList.push_back(auxiliaryRead(child)); else
      if (tag=="materialref") materialPtr = getMaterial(GenerateName(refRead(child))); else
      if (tag=="solidref") solidPtr = getSolid(GenerateName(refRead(child)));
   }

   pAssembly = 0;
   pMotherLogical = new G4LogicalVolume(solidPtr,materialPtr,GenerateName(name),0,0,0);

   if (!auxList.empty()) auxMap[pMotherLogical] = auxList;

   volume_contentRead(volumeElement);
}

void G4GDMLReadStructure::skinsurfaceRead(const xercesc::DOMElement* const skinsurfaceElement) {

   G4String name;
   G4LogicalVolume* logvol = 0;
   G4SurfaceProperty* prop = 0;

   const xercesc::DOMNamedNodeMap* const attributes = skinsurfaceElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="name") name = GenerateName(attValue); else
      if (attName=="surfaceproperty") prop = getSurfaceProperty(GenerateName(attValue));
   }

   for (xercesc::DOMNode* iter = skinsurfaceElement->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="volumeref") logvol = getVolume(GenerateName(refRead(child))); else
      G4Exception("GDML Reader: ERROR! Unknown tag in skinsurface: "+tag);
   }

   new G4LogicalSkinSurface(name,logvol,prop);
}

void G4GDMLReadStructure::volume_contentRead(const xercesc::DOMElement* const volumeElement) {

   for (xercesc::DOMNode* iter = volumeElement->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="auxiliary" || tag=="materialref" || tag=="solidref") { } else // These are already processed in volmeRead
      if (tag=="paramvol") paramvolRead(child,pMotherLogical); else
      if (tag=="physvol") physvolRead(child); else
      if (tag=="replicavol") replicavolRead(child); else
      if (tag=="divisionvol") divisionvolRead(child); else
      if (tag=="loop") loopRead(child,&G4GDMLRead::volume_contentRead); else
      G4Exception("GDML Reader: ERROR! Unknown tag in volume: "+tag);
   }
}

void G4GDMLReadStructure::structureRead(const xercesc::DOMElement* const structureElement) {

   for (xercesc::DOMNode* iter = structureElement->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="assembly") assemblyRead(child); else
      if (tag=="bordersurface") bordersurfaceRead(child); else
      if (tag=="skinsurface") skinsurfaceRead(child); else
      if (tag=="volume") volumeRead(child); else
      if (tag=="loop") loopRead(child,&G4GDMLRead::structureRead); else
      G4Exception("GDML Reader: ERROR! Unknown tag in structure: "+tag);
   }
}

G4AssemblyVolume* G4GDMLReadStructure::getAssembly(const G4String& ref) {

   if (assemblyMap.find(ref) == assemblyMap.end()) return 0;  // No error message is displayed!

   return assemblyMap[ref];
}

G4VPhysicalVolume* G4GDMLReadStructure::getPhysvol(const G4String& ref) const {

   G4VPhysicalVolume* physvolPtr = G4PhysicalVolumeStore::GetInstance()->GetVolume(ref,false);

   if (!physvolPtr) G4Exception("GDML Reader: ERROR! Referenced physvol '"+ref+"' was not found!");

   return physvolPtr;
}

G4LogicalVolume* G4GDMLReadStructure::getVolume(const G4String& ref) const {

   G4LogicalVolume *volumePtr = G4LogicalVolumeStore::GetInstance()->GetVolume(ref,false);

   if (!volumePtr) G4Exception("GDML Reader: ERROR! Referenced volume '"+ref+"' was not found!");

   return volumePtr;
}

G4GDMLAuxListType G4GDMLReadStructure::getVolumeAuxiliaryInformation(G4LogicalVolume* logvol) {

   if (auxMap.find(logvol) != auxMap.end()) return auxMap[logvol];
   else return G4GDMLAuxListType();
}
