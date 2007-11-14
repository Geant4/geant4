#include "G4GDMLStructure.hh"

G4GDMLStructure::G4GDMLStructure() {

   evaluator = G4GDMLEvaluator::GetInstance();

   parser = NULL;
}

bool G4GDMLStructure::directionRead(const xercesc::DOMElement* const element,EAxis& axis) {

   G4String x;
   G4String y;
   G4String z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="x") { x = attribute_value; }
      if (attribute_name=="y") { y = attribute_value; }
      if (attribute_name=="z") { z = attribute_value; }
   }

   G4double _x;
   G4double _y;
   G4double _z;

   if (!evaluator->Evaluate(_x,x)) return false;
   if (!evaluator->Evaluate(_y,y)) return false;
   if (!evaluator->Evaluate(_z,z)) return false;

   if (_x == 1.0 && _y == 0.0 && _z == 0.0) { axis = kXAxis; return true; } else
   if (_x == 0.0 && _y == 1.0 && _z == 0.0) { axis = kYAxis; return true; } else
   if (_x == 0.0 && _y == 0.0 && _z == 1.0) { axis = kZAxis; return true; }

   G4cout << "GDML ERROR! Only directions along axes are supported!"  << G4endl;

   return false;
}

bool G4GDMLStructure::divisionvolRead(const xercesc::DOMElement* const element,G4LogicalVolume* pMother) {

   G4String unit;
   G4String axis;
   G4String width;
   G4String offset;
   G4String number;
   G4String volumeref;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="unit"  ) { unit   = attribute_value; }
      if (attribute_name=="axis"  ) { unit   = attribute_value; }
      if (attribute_name=="width" ) { width  = attribute_value; }
      if (attribute_name=="offset") { offset = attribute_value; }
      if (attribute_name=="number") { number = attribute_value; }
   }

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="volumeref") { if (!refRead(child,volumeref)) return false; } 
   }

   G4LogicalVolume* pLogical = G4LogicalVolumeStore::GetInstance()->GetVolume(volumeref,false);

   if (pLogical == 0) {
   
      G4cout << "GDML ERROR! Referenced volume '" << volumeref << "' was not found in divisionvol in volume '" << pMother->GetName() << "'!" << G4endl;
      return false;
   }

   EAxis    _axis   = kZAxis;
   G4double _width  = 0.0;
   G4double _offset = 0.0;
   G4double _number = 0.0;
   
   if (axis=="kXAxis") { _axis = kXAxis; } else
   if (axis=="kYAxis") { _axis = kYAxis; } else
   if (axis=="kZAxis") { _axis = kZAxis; } else
   if (axis=="kRho"  ) { _axis = kRho;   } else
   if (axis=="kPhi"  ) { _axis = kPhi;   }

   if (!evaluator->Evaluate(_width ,width ,unit)) return false;
   if (!evaluator->Evaluate(_offset,offset,unit)) return false;
   if (!evaluator->Evaluate(_number,number))      return false;

   new G4PVDivision("",pLogical,pMother,_axis,(G4int)_number,_width,_offset);

   return true;
}

bool G4GDMLStructure::fileRead(const xercesc::DOMElement* const element,G4String& volname) {

   G4String name;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name"   ) { name    = attribute_value; } else
      if (attribute_name=="volname") { volname = attribute_value; }
   }

   G4String temp = module;  // push

   if (!gdmlRead(name,parser)) return false;

   if (!volname.empty()) return true;

   volname = setup.GetS("Default");

   volname = module + volname;

   module = temp; // pop

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

      if (tag=="file"       ) { if (!fileRead(child,volumeref  )) return false; } else
      if (tag=="volumeref"  ) { if (!refRead (child,volumeref  )) return false; } else
      if (tag=="positionref") { if (!refRead (child,positionref)) return false; } else
      if (tag=="rotationref") { if (!refRead (child,rotationref)) return false; }
   }

   volumeref = module + volumeref;

   G4LogicalVolume *pCurrentLogical = G4LogicalVolumeStore::GetInstance()->GetVolume(volumeref,false);

   if (pCurrentLogical == 0) {
   
      G4cout << "GDML: Error! Referenced volume '" << volumeref << "' was not found!" << G4endl;
      G4cout << "GDML: Error in physvol in volume '" << pMotherLogical->GetName() << "'!" << G4endl;
      return false;
   }

   G4ThreeVector tlate;

   if (!positionref.empty()) {
   
      G4ThreeVector *posPtr = solids.define.GetPosition(positionref);

      if (posPtr == 0) {

         G4cout << "GDML: Error in physvol in volume '" << pMotherLogical->GetName() << "'!" << G4endl;
         return false;
      }

      tlate = *posPtr;            
   }

   G4RotationMatrix *pRot=0;

   if (!rotationref.empty()) {

      G4ThreeVector *anglePtr = solids.define.GetRotation(rotationref);

      if (anglePtr == 0) {

         G4cout << "GDML: Error in physvol in volume '" << pMotherLogical->GetName() << "'!" << G4endl;
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

bool G4GDMLStructure::quantityRead(const xercesc::DOMElement* const element,G4double& _value) {

   G4String value;
   G4String unit;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="value") { value = attribute_value; }
      if (attribute_name=="unit" ) { unit  = attribute_value; }
   }

   return evaluator->Evaluate(_value,value,unit);
}

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

bool G4GDMLStructure::replicate_along_axisRead(const xercesc::DOMElement* const element,G4double& _width,G4double& _offset,EAxis& _axis) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="width"    ) { if (!quantityRead(child,_width))  return false; } else
      if (tag=="offset"   ) { if (!quantityRead(child,_offset)) return false; } else
      if (tag=="direction") { if (!directionRead(child,_axis))  return false; }
   }

   return true;
}

bool G4GDMLStructure::replicavolRead(const xercesc::DOMElement* const element,G4LogicalVolume* pMother) {

   G4String volumeref;
   G4String numb;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="numb") { numb = attribute_value; }
   }

   G4double _numb;
   G4double _width = 0.0;
   G4double _offset = 0.0;
   EAxis _axis;

   if (!evaluator->Evaluate(_numb,numb)) return false;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="volumeref"           ) { if (!refRead(child,volumeref))                             return false; } else
      if (tag=="replicate_along_axis") { if (!replicate_along_axisRead(child,_width,_offset,_axis)) return false; }
   }

   G4LogicalVolume* pLogical = G4LogicalVolumeStore::GetInstance()->GetVolume(volumeref,false);

   if (pLogical == 0) {
   
      G4cout << "GDML ERROR! Referenced volume '" << volumeref << "' was not found in replicavol in volume '" << pMother->GetName() << "'!" << G4endl;
      return false;
   }

   new G4PVReplica("",pLogical,pMother,_axis,(G4int)_numb,_width,_offset);

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

   G4Material *materialPtr = materials.Get(materialref);
   G4VSolid* solidPtr = solids.Get(solidref);

   if (!solidPtr || !materialPtr) {
   
      G4cout << "GDML: Error in volume '" << name << "'!" << G4endl;
      return false;
   }
 
   G4LogicalVolume *volumePtr = new G4LogicalVolume(solidPtr,materialPtr,module+name,0,0,0);
   
   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
  
      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="physvol"    ) { if (!physvolRead    (child,volumePtr)) return false; } 
      if (tag=="replicavol" ) { if (!replicavolRead (child,volumePtr)) return false; }
      if (tag=="divisionvol") { if (!divisionvolRead(child,volumePtr)) return false; }
   }

   return true;
}

bool G4GDMLStructure::Read(const xercesc::DOMElement* const element,const G4String& newModule) {

   module = newModule;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="volume") { if (!volumeRead(child)) return false; } else
      {
	 G4cout << "GDML ERROR! Unknown tag in structure: " << tag << G4endl;
         return false;
      }
   }
  
   return true;
}

bool G4GDMLStructure::gdmlRead(const G4String& fileName,xercesc::XercesDOMParser* newParser) {

   G4String newModule;

   if (parser != NULL) newModule = fileName + "::"; // This must be an external module

   parser = newParser;

   try {

      parser->parse(fileName.c_str());
   }
   catch (const xercesc::XMLException &e) {
   
      char* message = xercesc::XMLString::transcode(e.getMessage());
      G4cout << "XML: " << message << G4endl;
      xercesc::XMLString::release(&message);
   }
   catch (const xercesc::DOMException &e) {
   
      char* message = xercesc::XMLString::transcode(e.getMessage());
      G4cout << "DOM: " << message << G4endl;
      xercesc::XMLString::release(&message);
   }
   catch (...) {

      G4cout << "Unexpected exception!" << G4endl;
   }

   xercesc::DOMDocument* doc = parser->getDocument();

   if (!doc) {
   
      G4cout << "Unable to open document '" << fileName << "'!" << G4endl;
      return false;
   }

   xercesc::DOMElement* element = doc->getDocumentElement();

   if (!element) {
   
      G4cout << "Empty document!" << G4endl;
      return false;
   }

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="define"   ) { if (!solids.define.Read(child,newModule)) return false; } else
      if (tag=="materials") { if (!materials.Read    (child,newModule)) return false; } else
      if (tag=="solids"   ) { if (!solids.Read       (child,newModule)) return false; } else
      if (tag=="setup"    ) { if (!setup.Read        (child,newModule)) return false; } else
      if (tag=="structure") { if (!Read              (child,newModule)) return false; }
   }

   return true;
}
