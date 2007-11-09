#include "G4GDMLStructure.hh"

bool G4GDMLStructure::refRead(const xercesc::DOMElement* const element,G4String& ref) {

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="ref") { ref = attribute_value; }
   }

   return true;
}

bool G4GDMLStructure::physvolRead(const xercesc::DOMElement* const element,G4LogicalVolume *pMotherLogical) {

   G4String volumeref;
   G4String positionref;
   G4String rotationref;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="volumeref"  ) { if (!refRead(child,volumeref  )) return false; } else
      if (tag=="positionref") { if (!refRead(child,positionref)) return false; } else
      if (tag=="rotationref") { if (!refRead(child,rotationref)) return false; } else
      {
	 G4cout << "GDML ERROR! Unsupported tag in physvol: " << tag << G4endl;
         return false;
      }
   }

   G4LogicalVolume *pCurrentLogical = G4LogicalVolumeStore::GetInstance()->GetVolume(volumeref,false);

   if (pCurrentLogical == 0) {
   
      G4cout << "G4GDMLParser: Error in the physvol of volume '" << pMotherLogical->GetName() << "': ";
      G4cout << "referenced volume '" << volumeref << "' can not be found!" << G4endl;
      return false;
   }

   G4ThreeVector tlate;

   if (positionref != "") {
   
      G4ThreeVector *posPtr = solids.define.GetPosition(positionref);

      if (posPtr == 0) {
      
         G4cout << "G4GDMLParser: Error in the physvol of volume '" << pMotherLogical->GetName() << "': ";
         G4cout << "referenced position '" << positionref << "' can not be found!" << G4endl;
         return false;
      }

      tlate = *posPtr;            
   }

   G4RotationMatrix *pRot=0;

   if (rotationref != "") {

      G4ThreeVector *anglePtr = solids.define.GetRotation(rotationref);

      if (anglePtr == 0) {
      
         G4cout << "G4GDMLParser: Error in the physvol of volume '" << pMotherLogical->GetName() << "': ";
         G4cout << "referenced rotation '" << rotationref << "' can not be found!" << G4endl;
         return false;
      }

      pRot = new G4RotationMatrix();

      pRot->rotateX(anglePtr->x());
      pRot->rotateY(anglePtr->y());
      pRot->rotateZ(anglePtr->z());
   }
 
   new G4PVPlacement(pRot,tlate,pCurrentLogical,"",pMotherLogical,false,0);
 
   return true;
}

bool G4GDMLStructure::volumeRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String solidref;
   G4String materialref;

   XMLCh *name_attr = xercesc::XMLString::transcode("name");
   name = xercesc::XMLString::transcode(element->getAttribute(name_attr));
   xercesc::XMLString::release(&name_attr);

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="materialref") { if (!refRead(child,materialref)) return false; } else
      if (tag=="solidref"   ) { if (!refRead(child,solidref   )) return false; }
   }

   G4VSolid *solidPtr = G4SolidStore::GetInstance()->GetSolid(solidref,false); 

   if (solidPtr == 0) {
   
      G4cout << "G4GDMLParser: Error in volume '" << name << "': ";
      G4cout << "referenced solid '" << solidref << "' can not be found!" << G4endl;
      return false;
   }
 
   G4Material *materialPtr = G4Material::GetMaterial(materialref,false);

   if (materialPtr == 0) {
   
      G4cout << "G4GDMLParser: Error in volume '" << name << "': ";
      G4cout << "referenced material '" << materialref << "' can not be found!" << G4endl;
      return false;
   }

   G4LogicalVolume *volumePtr = new G4LogicalVolume(solidPtr,materialPtr,name,0,0,0); // Create the logical volume, so that the physvol children can be appended!

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
  
      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="physvol") { if (!physvolRead(child,volumePtr)) return false; } 
   }

   return true;
}

bool G4GDMLStructure::Read(const xercesc::DOMElement* const element) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="volume") { if (!volumeRead(child)) return false; } else
      {
	 G4cout << "GDML ERROR! Unsupported tag in structure: " << tag << G4endl;
         return false;
      }
   }
  
   return true;
}
