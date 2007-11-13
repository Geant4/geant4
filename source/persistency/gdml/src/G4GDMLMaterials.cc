#include "G4GDMLMaterials.hh"

G4GDMLMaterials::G4GDMLMaterials() {

   evaluator = G4GDMLEvaluator::GetInstance();
}

bool G4GDMLMaterials::atomRead(const xercesc::DOMElement* const element,double& _value) {

   std::string value;
   std::string unit = "g/mole";

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="value") { value = attribute_value; } else
      if (attribute_name=="unit ") { unit  = attribute_value; } else
      {
      }
   }

   return evaluator->Evaluate(_value,value,unit);
}

bool G4GDMLMaterials::DRead(const xercesc::DOMElement* const element,double& _value) {

   std::string value;
   std::string unit = "g/cm3";

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="value") { value = attribute_value; } else
      if (attribute_name=="unit ") { unit  = attribute_value; } else
      {
      }
   }

   return evaluator->Evaluate(_value,value,unit);
}

bool G4GDMLMaterials::elementRead(const xercesc::DOMElement* const element) {

   std::string name;
   std::string formula;
   std::string Z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name"   ) { name    = attribute_value; } else
      if (attribute_name=="formula") { formula = attribute_value; } else
      if (attribute_name=="Z"      ) { Z       = attribute_value; } else
      {
      }
   }

   double _Z;
   double _a;

   if (!evaluator->Evaluate(_Z,Z)) return false;

   int nComponents = 0;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="atom"    ) { if (!atomRead(child,_a)) return false; } else
      if (tag=="fraction") { nComponents++;                         } else
      {
         std::cout << "ERROR! Unsupported tag in element: " << tag << std::endl;
         return false;
      }
   }

   if (nComponents > 0) return mixtureRead(element,new G4Element(module+name,formula,nComponents));

   new G4Element(module+name,formula,_Z,_a);

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

bool G4GDMLMaterials::isotopeRead(const xercesc::DOMElement* const element) {

   std::string name;
   std::string Z;
   std::string N;

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
      if (attribute_name=="N"   ) { N    = attribute_value; } else
      {
      }
   }

   double _Z;
   double _N;
   double _a;

   if (!evaluator->Evaluate(_Z,Z)) return false;
   if (!evaluator->Evaluate(_N,N)) return false;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="atom") { if (!atomRead(child,_a)) return false; }
   }

   new G4Isotope(module+name,(G4int)_Z,(G4int)_N,_a);

   return true;
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
      if (attribute_name=="Z"   ) { Z    = attribute_value; }
   }
  
   double _Z;
   double _D;
   double _a;

   if (!evaluator->Evaluate(_Z,Z)) return false;

   int nComponents = 0;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="D"       ) { if (!DRead   (child,_D)) return false; } else
      if (tag=="atom"    ) { if (!atomRead(child,_a)) return false; } else
      if (tag=="fraction") { nComponents++;                         } else
      {
	 G4cout << "GDML ERROR! Unsupported tag in material: " << tag << G4endl;
         return false;
      }
   }

   if (nComponents > 0) return mixtureRead(element,new G4Material(module+name,_D,nComponents));

   new G4Material(module+name,_Z,_a,_D);

   return true;
}

bool G4GDMLMaterials::mixtureRead(const xercesc::DOMElement *const element,G4Element *ele) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="fraction") {

         std::string ref;
	 double _n;
	 
	 if (!fractionRead(child,_n,ref)) return false;

         G4Isotope *isotopePtr = G4Isotope::GetIsotope(ref,false);

         if (isotopePtr == 0) {
   
            G4cout << "GDML ERROR! Referenced isotope '" << ref << "' in element '" << ele->GetName() << "' was not found!" << G4endl;   
            return false;
         }      
      
         ele->AddIsotope(isotopePtr,_n);
      } else {
      
         std::cout << "ERROR! Unsupported tag in mixture element: " << tag << std::endl;
         return false;
      }
   }

   return true;
}

bool G4GDMLMaterials::mixtureRead(const xercesc::DOMElement *const element,G4Material *material) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="D"       ) { /*already processed*/ } else
      if (tag=="fraction") {

         std::string ref;
	 double _n;
	 
	 if (!fractionRead(child,_n,ref)) return false;

         G4Material *materialPtr = G4Material::GetMaterial(ref,false);
         G4Element *elementPtr = G4Element::GetElement(ref,false);

         if (materialPtr != 0) material->AddMaterial(materialPtr,_n); else
	 if (elementPtr != 0) material->AddElement(elementPtr,_n);

         if ((materialPtr == 0) && (elementPtr == 0)) {
   
            G4cout << "GDML ERROR! Referenced material/element '" << ref << "' in material '" << material->GetName() << "' was not found!" << G4endl;   
            return false;
         }      
      } else {
      
         std::cout << "ERROR! Unsupported tag in mixture material: " << tag << std::endl;
         return false;
      }
   }

   return true;
}

bool G4GDMLMaterials::Read(const xercesc::DOMElement* const element,const G4String& newModule) {

   module = newModule;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != NULL;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="element" ) { if (!elementRead (child)) return false; } else 
      if (tag=="isotope" ) { if (!isotopeRead (child)) return false; } else 
      if (tag=="material") { if (!materialRead(child)) return false; } else 
      {
	 G4cout << "GDML ERROR! Unsupported tag in materials: " << tag << G4endl;
         return false;
      }
   }

   return true;
}
