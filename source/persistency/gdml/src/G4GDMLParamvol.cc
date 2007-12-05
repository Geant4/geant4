#include "G4GDMLParamvol.hh"

void G4GDMLParamvol::box_parametersRead(const xercesc::DOMElement* const element,G4GDMLParameterisation::PARAMETER& parameter) {

   G4String lunit("1");
   G4String x;
   G4String y;
   G4String z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="lunit") lunit = attribute_value; else
      if (attribute_name=="x") x = attribute_value; else
      if (attribute_name=="y") y = attribute_value; else
      if (attribute_name=="z") z = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);

   parameter.dimension[0] = eval.Evaluate(x)*_lunit;
   parameter.dimension[1] = eval.Evaluate(y)*_lunit;
   parameter.dimension[2] = eval.Evaluate(z)*_lunit;

   parameter.dimension[0] *= 0.5;
   parameter.dimension[1] *= 0.5;
   parameter.dimension[2] *= 0.5;
}

G4String G4GDMLParamvol::refRead(const xercesc::DOMElement* const element) {

   G4String ref;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="ref") ref = attribute_value;
   }

   return ref;
}

void G4GDMLParamvol::parametersRead(const xercesc::DOMElement* const element) {

   G4ThreeVector rotation;
   G4ThreeVector position;

   G4GDMLParameterisation::PARAMETER parameter;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());
      
      if (tag=="rotation") rotation = rotationRead(child); else
      if (tag=="position") position = positionRead(child); else
      if (tag=="box_parameters") box_parametersRead(child,parameter);
   }

   parameter.pRot = new G4RotationMatrix();
   
   parameter.pRot->rotateX(rotation.x());
   parameter.pRot->rotateY(rotation.y());
   parameter.pRot->rotateZ(rotation.z());

   parameter.position = position;

   parameterisation->addParameter(parameter);
}

void G4GDMLParamvol::paramvolRead(const xercesc::DOMElement* const element,G4LogicalVolume* mother) {

   G4String volumeref;

   parameterisation = new G4GDMLParameterisation();

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

       if (tag=="parameters") parametersRead(child);
       if (tag=="volumeref") volumeref = refRead(child);
    }

   G4LogicalVolume* logvol = getVolume(volumeref);

   new G4PVParameterised("",logvol,mother,kZAxis,parameterisation->getSize(),parameterisation);
}
