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

G4GDMLReadStructure::AuxPairType G4GDMLReadStructure::auxiliaryRead(const xercesc::DOMElement* const element) {

   G4String auxtype;
   G4String auxvalue;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="auxtype") auxtype = attValue; else
      if (attName=="auxvalue") auxvalue = attValue;
   }

   return AuxPairType(auxtype,auxvalue);
}

void G4GDMLReadStructure::divisionvolRead(const xercesc::DOMElement* const divisionvolElement) {

   G4double unit = 1.0;
   G4double width = 0.0;
   G4double offset = 0.0;
   G4int number = 0;
   G4String volumeref;
   EAxis axis = kUndefined;

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

      if (tag=="volumeref") volumeref = refRead(child);
   }

   G4LogicalVolume* pLogical = getVolume(GenerateName(volumeref));

   G4String name = pLogical->GetName() + "_in_" + pMotherLogical->GetName();

   if (number != 0 && width == 0.0) new G4PVDivision(name,pLogical,pMotherLogical,axis,number,offset); else
   if (number == 0 && width != 0.0) new G4PVDivision(name,pLogical,pMotherLogical,axis,width,offset); else
   if (number != 0 && width != 0.0) new G4PVDivision(name,pLogical,pMotherLogical,axis,number,width,offset); else
   G4Exception("GDML Reader: ERROR! Both 'number' and 'width' are zeros in divisionvol: "+name);
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
      if (tag=="position") position = vectorRead(child); else
      if (tag=="rotation") rotation = vectorRead(child); else
      if (tag=="scale") scale = vectorRead(child); else
      if (tag=="positionref") position = *getPosition(GenerateName(refRead(child))); else
      if (tag=="rotationref") rotation = *getRotation(GenerateName(refRead(child))); else
      if (tag=="scaleref") scale = *getScale(GenerateName(refRead(child))); else
      G4Exception("GDML Reader: ERROR! Unknown tag in physvol: "+tag);
   }

   G4Transform3D transform(getRotationMatrix(rotation).inverse(),position);
   transform = transform*G4Scale3D(scale.x(),scale.y(),scale.z());

   G4PhysicalVolumesPair pvPair = G4ReflectionFactory::Instance()->Place(transform,name,logvol,pMotherLogical,false,0,false);

   if (name.empty()) {

      if (pvPair.first != 0) GeneratePhysvolName(pvPair.first);
      if (pvPair.second != 0) GeneratePhysvolName(pvPair.second);
   }
}

void G4GDMLReadStructure::replicavolRead(const xercesc::DOMElement* const replicavolElement) {

   G4double unit = 1.0;
   G4double width = 0.0;
   G4double offset = 0.0;
   G4int number = 0;
   G4String volumeref;
   EAxis axis = kUndefined;

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

      if (tag=="volumeref") volumeref = refRead(child);
   }

   G4LogicalVolume* pLogical = getVolume(GenerateName(volumeref));

   G4String name = pLogical->GetName() + "_in_" + pMotherLogical->GetName();

   new G4PVReplica(name,pLogical,pMotherLogical,axis,number,width,offset);
}

void G4GDMLReadStructure::volumeRead(const xercesc::DOMElement* const volumeElement) {

   G4String name;

   G4VSolid* solidPtr = 0;
   G4Material* materialPtr = 0;

   AuxListType auxList;

   XMLCh *name_attr = xercesc::XMLString::transcode("name");
   name = xercesc::XMLString::transcode(volumeElement->getAttribute(name_attr));
   xercesc::XMLString::release(&name_attr);

   for (xercesc::DOMNode* iter = volumeElement->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="auxiliary") auxList.push_back(auxiliaryRead(child)); else
      if (tag=="materialref") materialPtr = getMaterial(GenerateName(refRead(child))); else
      if (tag=="solidref") solidPtr = getSolid(GenerateName(refRead(child)));
   }

   pMotherLogical = new G4LogicalVolume(solidPtr,materialPtr,GenerateName(name),0,0,0);

   if (!auxList.empty()) auxMap[pMotherLogical] = auxList;

   const G4LogicalVolumeStore* volumeList = G4LogicalVolumeStore::GetInstance();   
   const size_t volumeCount = volumeList->size();

   volume_contentRead(volumeElement);

   if (volumeCount != volumeList->size()) {   // New logical volume can come from "volume_contentRead()".

      volumeList->DeRegister(pMotherLogical);  // "pMotherLogical" must be the last in the list, since the new volume can be referenced!
      volumeList->Register(pMotherLogical);    // This arranging only matters if multiple GDML files are written out into a single GDML file.
   }
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

      if (tag=="volume") volumeRead(child); else      
      if (tag=="loop") loopRead(child,&G4GDMLRead::structureRead); else
      G4Exception("GDML Reader: ERROR! Unknown tag in structure: "+tag);
   }
}

G4LogicalVolume* G4GDMLReadStructure::getVolume(const G4String& ref) const {

   G4LogicalVolume *volumePtr = G4LogicalVolumeStore::GetInstance()->GetVolume(ref,false);

   if (!volumePtr) G4Exception("GDML Reader: ERROR! Referenced volume '"+ref+"' was not found!");

   return volumePtr;
}

G4GDMLReadStructure::AuxListType G4GDMLReadStructure::getVolumeAuxiliaryInformation(const G4LogicalVolume* const ptr) {

     if (auxMap.find(ptr) != auxMap.end()) return auxMap[ptr];
     else return G4GDMLReadStructure::AuxListType();
}

const G4GDMLReadStructure::AuxMapType* G4GDMLReadStructure::getAuxiliaryMap() {

   return &auxMap;
}
