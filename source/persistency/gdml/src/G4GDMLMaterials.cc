#include "G4GDMLMaterials.hh"

G4GDMLMaterials::G4GDMLMaterials() {

   evaluator = G4GDMLEvaluator::GetInstance();
}

bool G4GDMLMaterials::compositeRead(const xercesc::DOMElement *const element,G4Material *material) {

   const xercesc::DOMNodeList* const children = element->getChildNodes();
   const  XMLSize_t nodeCount = children->getLength();

   for (XMLSize_t i=0;i<nodeCount;i++) {

      xercesc::DOMNode* node = children->item(i);
      
      if (node->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(node);   

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="fraction") {

         std::string ref;
	 double _n;
	 
	 if (!fractionRead(child,_n,ref)) return false;

         G4Material *refPtr = G4Material::GetMaterial(ref,false);

         if (refPtr == 0) {
   
            std::cout << std::endl;
            std::cout << "GDML ERROR! Referenced material '" << ref << "' in material '" << material->GetName() << "' was not found!" << std::endl;   
            std::cout << std::endl;
            return false;
         }

         material->AddMaterial(refPtr,_n);            
      }
   }

   return true;
}

bool G4GDMLMaterials::elementRead(const xercesc::DOMElement* const) {

   return true;
}

bool G4GDMLMaterials::fractionRead(const xercesc::DOMElement* const element,double& _n,std::string& ref) {

   std::string n;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="n"  ) { n   = attribute_value; } else
      if (attribute_name=="ref") { ref = attribute_value; } else
      {
      }
   }

   return evaluator->Evaluate(_n,n);
}

bool G4GDMLMaterials::materialRead(const xercesc::DOMElement* const element) {

   std::string name;
   std::string Z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") { name = attribute_value; } else
      if (attribute_name=="Z"   ) { Z    = attribute_value; } else
      {
      }
   }
  
   double _Z;
   double _D;
   double _a;

   if (!evaluator->Evaluate(_Z,Z)) return false;

   int nComponents = 0;

   const xercesc::DOMNodeList* const children = element->getChildNodes();
   XMLSize_t nodeCount = children->getLength();

   for (XMLSize_t i=0;i<nodeCount;i++) {

      xercesc::DOMNode* node = children->item(i);
      
      if (node->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(node);  

      const std::string material_tag = xercesc::XMLString::transcode(child->getTagName());

      if (material_tag=="D"       ) { if (!valueRead(child,_D)) return false; } else
      if (material_tag=="atom"    ) { if (!valueRead(child,_a)) return false; } else
      if (material_tag=="fraction") { nComponents++;                          } else
      if (material_tag=="element" ) { nComponents++;                          } else
      {
      }
   }

   _D *= g/cm3;
   _a *= g/mole;

   if (nComponents > 0)
      return compositeRead(element,new G4Material(name,_D,nComponents));

   new G4Material(name,_Z,_a,_D);

   return true;
}

bool G4GDMLMaterials::valueRead(const xercesc::DOMElement* const element,double& _value) {

   std::string value;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="value") { value = attribute_value; } else
      {
      }
   }

   return evaluator->Evaluate(_value,value);
}

bool G4GDMLMaterials::Read(const xercesc::DOMElement* const element) {

   const xercesc::DOMNodeList* const children = element->getChildNodes();
   XMLSize_t elementCount = children->getLength();

   for (XMLSize_t element_index=0;element_index<elementCount;element_index++) {

      xercesc::DOMNode* element_node = children->item(element_index);
      
      if (element_node->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;
   
      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(element_node);   

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="element" ) { if (!elementRead (child)) return false; } else 
      if (tag=="material") { if (!materialRead(child)) return false; } else 
      {
         std::cout << std::endl;
	 std::cout << "GDML ERROR! Unsupported tag in materials: " << tag << std::endl;
         std::cout << std::endl;
         return false;
      }
   }

   return true;
}
