#include "G4GDMLDefine.hh"

G4GDMLDefine::G4GDMLDefine() {

    evaluator = G4GDMLEvaluator::GetInstance();
}

G4GDMLDefine::~G4GDMLDefine() {
}

bool G4GDMLDefine::constantRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String value;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name" ) { name  = attribute_value; } else
      if (attribute_name=="value") { value = attribute_value; }
   }

   G4double _value;

   if (!evaluator->Evaluate(_value,value)) return false;

   return evaluator->RegisterConstant(name,_value);
}

bool G4GDMLDefine::positionRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String unit;
   G4String x;
   G4String y;
   G4String z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") { name = attribute_value; } else
      if (attribute_name=="unit") { unit = attribute_value; } else
      if (attribute_name=="x"   ) { x    = attribute_value; } else
      if (attribute_name=="y"   ) { y    = attribute_value; } else
      if (attribute_name=="z"   ) { z    = attribute_value; }
   }

   G4double _x;
   G4double _y;
   G4double _z;

   if (!evaluator->Evaluate(_x,x,unit)) return false;
   if (!evaluator->Evaluate(_y,y,unit)) return false;
   if (!evaluator->Evaluate(_z,z,unit)) return false;

   positionMap[module+name] = new G4ThreeVector(_x,_y,_z);

   return true;
}

bool G4GDMLDefine::rotationRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String unit;
   G4String x;
   G4String y;
   G4String z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") { name = attribute_value; } else
      if (attribute_name=="unit") { unit = attribute_value; } else
      if (attribute_name=="x"   ) { x    = attribute_value; } else
      if (attribute_name=="y"   ) { y    = attribute_value; } else
      if (attribute_name=="z"   ) { z    = attribute_value; }
   }
   
   G4double _x;
   G4double _y;
   G4double _z;

   if (!evaluator->Evaluate(_x,x,unit)) return false;
   if (!evaluator->Evaluate(_y,y,unit)) return false;
   if (!evaluator->Evaluate(_z,z,unit)) return false;

   rotationMap[module+name] = new G4ThreeVector(_x,_y,_z);

   return true;
}

bool G4GDMLDefine::Read(const xercesc::DOMElement* const element,const G4String& newModule) {

   module = newModule;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="constant") { if (!constantRead(child)) return false; } else
      if (tag=="position") { if (!positionRead(child)) return false; } else
      if (tag=="rotation") { if (!rotationRead(child)) return false; } else
      {
	 G4cout << "GDML: Error! Unknown tag in define: " << tag << G4endl;
         return false;
      }
   }
            
   return true;
}

G4ThreeVector* G4GDMLDefine::GetPosition(const G4String& ref) {

   G4String full_ref = module + ref;

   if (positionMap.find(full_ref) != positionMap.end()) return positionMap[full_ref];

   G4cout << "GDML: Error! Referenced position '" << full_ref << "' was not found!" << G4endl;

   return 0;
}

G4ThreeVector* G4GDMLDefine::GetRotation(const G4String& ref) {

   G4String full_ref = module + ref;

   if (rotationMap.find(full_ref) != rotationMap.end()) return rotationMap[full_ref];

   G4cout << "GDML: Error! Referenced rotation '" << full_ref << "' was not found!" << G4endl;

   return 0;
}
