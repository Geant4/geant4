#include "G4GDMLStructure.hh"

bool G4GDMLStructure::physvolRead(const xercesc::DOMElement* const element,G4LogicalVolume *pMotherLogical) {

   std::string volumeref;
   std::string positionref;
   std::string rotationref;
   std::string scaleref;

   XMLCh *ref_attr = xercesc::XMLString::transcode("ref");

   const xercesc::DOMNodeList* const children = element->getChildNodes();
   const  XMLSize_t nodeCount = children->getLength();

   for (XMLSize_t i=0;i<nodeCount;i++) {

      xercesc::DOMNode* node = children->item(i);
      
      if (node->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;
   
      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(node);   

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="volumeref"  ) { volumeref   = xercesc::XMLString::transcode(child->getAttribute(ref_attr)); } else
      if (tag=="positionref") { positionref = xercesc::XMLString::transcode(child->getAttribute(ref_attr)); } else
      if (tag=="rotationref") { rotationref = xercesc::XMLString::transcode(child->getAttribute(ref_attr)); } else
      if (tag=="scaleref"   ) { scaleref    = xercesc::XMLString::transcode(child->getAttribute(ref_attr)); }
   }

   xercesc::XMLString::release(&ref_attr);

   G4LogicalVolume *pCurrentLogical = G4LogicalVolumeStore::GetInstance()->GetVolume(volumeref,false);

   if (pCurrentLogical == 0) {
   
      std::cout << "G4GDMLParser: Error in the physvol of volume '" << pMotherLogical->GetName() << "': ";
      std::cout << "referenced volume '" << volumeref << "' can not be found!" << std::endl;
      return false;
   }

   G4ThreeVector tlate;

   if (positionref != "") {
   
      G4ThreeVector *posPtr = solids.define.GetPosition(positionref);

      if (posPtr == 0) {
      
         std::cout << "G4GDMLParser: Error in the physvol of volume '" << pMotherLogical->GetName() << "': ";
         std::cout << "referenced position '" << positionref << "' can not be found!" << std::endl;
         return false;
      }

      tlate = *posPtr;            
   }

   G4RotationMatrix *pRot=0;

   if (rotationref != "") {

      G4ThreeVector *anglePtr = solids.define.GetRotation(rotationref);

      if (anglePtr == 0) {
      
         std::cout << "G4GDMLParser: Error in the physvol of volume '" << pMotherLogical->GetName() << "': ";
         std::cout << "referenced rotation '" << rotationref << "' can not be found!" << std::endl;
         return false;
      }

      pRot = new G4RotationMatrix();

      pRot->rotateX(anglePtr->x());
      pRot->rotateY(anglePtr->y());
      pRot->rotateZ(anglePtr->z());
   }
 
   new G4PVPlacement(pRot,tlate,pCurrentLogical,"",pMotherLogical,false,0);  // Name does not matter, let it be an empty string
 
   return true;
}

bool G4GDMLStructure::volumeRead(const xercesc::DOMElement* const element) {

   std::string name;
   std::string solidref;
   std::string materialref;

   XMLCh *name_attr = xercesc::XMLString::transcode("name");
   name = xercesc::XMLString::transcode(element->getAttribute(name_attr));
   xercesc::XMLString::release(&name_attr);

   const xercesc::DOMNodeList* const children = element->getChildNodes();
   XMLSize_t nodeCount = children->getLength();

   for (XMLSize_t i=0;i<nodeCount;i++) {

      xercesc::DOMNode* node = children->item(i);
      
      if (node->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;
   
      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(node);   
      
      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      XMLCh *ref_attr = xercesc::XMLString::transcode("ref");

      if (tag=="materialref") { materialref = xercesc::XMLString::transcode(child->getAttribute(ref_attr)); } else
      if (tag=="solidref"   ) { solidref    = xercesc::XMLString::transcode(child->getAttribute(ref_attr)); }

      xercesc::XMLString::release(&ref_attr);
   }

   G4VSolid *solidPtr = G4SolidStore::GetInstance()->GetSolid(solidref,false); 

   if (solidPtr == 0) {
   
      std::cout << "G4GDMLParser: Error in volume '" << name << "': ";
      std::cout << "referenced solid '" << solidref << "' can not be found!" << std::endl;
      return false;
   }
 
   G4Material *materialPtr = G4Material::GetMaterial(materialref,false);

   if (materialPtr == 0) {
   
      std::cout << "G4GDMLParser: Error in volume '" << name << "': ";
      std::cout << "referenced material '" << materialref << "' can not be found!" << std::endl;
      return false;
   }

   // At first, the logical volume must be created. The contained physical volumes can be created afterwards!

   G4LogicalVolume *volumePtr = new G4LogicalVolume(solidPtr,materialPtr,name,0,0,0);

   for (XMLSize_t i=0;i<nodeCount;i++) {

      xercesc::DOMNode* node = children->item(i);
      
      if (node->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;
   
      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(node);   
      
      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="physvol") { if (!physvolRead(child,volumePtr)) return false; } 
   }

   return true;
}

bool G4GDMLStructure::Read(const xercesc::DOMElement* const element) {

   const xercesc::DOMNodeList* const children = element->getChildNodes();
   const  XMLSize_t nodeCount = children->getLength();

   for (XMLSize_t i=0;i<nodeCount;i++) {

      xercesc::DOMNode* node = children->item(i);
      
      if (node->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;
   
      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(node);   

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="volume") { if (!volumeRead(child)) return false; } else
      {
	 std::cout << "   G4GDMLParser: Error! Tag '" << tag << "' is not recognized in structure!" << std::endl;
         return false;
      }
   }

   return true;
}
