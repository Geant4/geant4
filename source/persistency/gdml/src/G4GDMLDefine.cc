#include "G4GDMLDefine.hh"

G4GDMLDefine::G4GDMLDefine() {

    evaluator = G4GDMLEvaluator::GetInstance();
}

G4GDMLDefine::~G4GDMLDefine() {
}

bool G4GDMLDefine::constantRead(const xercesc::DOMElement* const element) {

   std::string name;
   std::string value;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name" ) { name  = attribute_value; } else
      if (attribute_name=="value") { value = attribute_value; } else
      {
         std::cout << std::endl;
         std::cout << "GDML ERROR! Unsupported attribute in constant '" << name << "': " << attribute_name << std::endl;
         std::cout << std::endl;
         return false;
      }
   }

   double _value;

   if (!evaluator->Evaluate(_value,value)) return false;

   return evaluator->RegisterConstant(name,_value);
}

bool G4GDMLDefine::positionRead(const xercesc::DOMElement* const element) {

   std::string unit;
   std::string name;
   std::string x;
   std::string y;
   std::string z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="unit") { unit = attribute_value; } else
      if (attribute_name=="name") { name = attribute_value; } else
      if (attribute_name=="x"   ) { x    = attribute_value; } else
      if (attribute_name=="y"   ) { y    = attribute_value; } else
      if (attribute_name=="z"   ) { z    = attribute_value; } else
      {
      }
   }

   double _x;
   double _y;
   double _z;

   if (!evaluator->Evaluate(_x,x,unit)) return false;
   if (!evaluator->Evaluate(_y,y,unit)) return false;
   if (!evaluator->Evaluate(_z,z,unit)) return false;

   positionMap[name] = new G4ThreeVector(_x,_y,_z);

   return true;
}

bool G4GDMLDefine::rotationRead(const xercesc::DOMElement* const element) {

   std::string unit;
   std::string name;
   std::string x;
   std::string y;
   std::string z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="unit") { unit = attribute_value; } else
      if (attribute_name=="name") { name = attribute_value; } else
      if (attribute_name=="x"   ) { x    = attribute_value; } else
      if (attribute_name=="y"   ) { y    = attribute_value; } else
      if (attribute_name=="z"   ) { z    = attribute_value; } else
      {
      }
   }
   
   double _x;
   double _y;
   double _z;

   if (!evaluator->Evaluate(_x,x,unit)) return false;
   if (!evaluator->Evaluate(_y,y,unit)) return false;
   if (!evaluator->Evaluate(_z,z,unit)) return false;

   rotationMap[name] = new G4ThreeVector(_x,_y,_z);

   return true;
}

bool G4GDMLDefine::Read(const xercesc::DOMElement* const element) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="constant") { if (!constantRead(child)) return false; } else
      if (tag=="position") { if (!positionRead(child)) return false; } else
      if (tag=="rotation") { if (!rotationRead(child)) return false; } else
      {
         std::cout << std::endl;
	 std::cout << "GDML ERROR! Unsupported tag in define: " << tag << std::endl;
         std::cout << std::endl;
         return false;
      }
   }
            
   return true;
}

G4ThreeVector *G4GDMLDefine::GetPosition(const std::string& ref) {

   return (positionMap.find(ref) == positionMap.end()) ? (0) : (positionMap[ref]);
}

G4ThreeVector *G4GDMLDefine::GetRotation(const std::string& ref) {

   return (rotationMap.find(ref) == rotationMap.end()) ? (0) : (rotationMap[ref]);
}
