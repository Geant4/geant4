#include "G4GDMLReadParamvol.hh"

void G4GDMLReadParamvol::box_dimensionsRead(const xercesc::DOMElement* const element,G4GDMLParameterisation::PARAMETER& parameter) {

   G4double lunit = 1.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="lunit") lunit = eval.Evaluate(attValue); else
      if (attName=="x") parameter.dimension[0] = eval.Evaluate(attValue); else
      if (attName=="y") parameter.dimension[1] = eval.Evaluate(attValue); else
      if (attName=="z") parameter.dimension[2] = eval.Evaluate(attValue);
   }

   parameter.dimension[0] *= 0.5*lunit;
   parameter.dimension[1] *= 0.5*lunit;
   parameter.dimension[2] *= 0.5*lunit;
}

void G4GDMLReadParamvol::trd_dimensionsRead(const xercesc::DOMElement* const element,G4GDMLParameterisation::PARAMETER& parameter) {

   G4double lunit = 1.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="lunit") lunit = eval.Evaluate(attValue); else
      if (attName=="x1") parameter.dimension[0] = eval.Evaluate(attValue); else
      if (attName=="x2") parameter.dimension[1] = eval.Evaluate(attValue); else
      if (attName=="y1") parameter.dimension[2] = eval.Evaluate(attValue); else
      if (attName=="y2") parameter.dimension[3] = eval.Evaluate(attValue); else
      if (attName=="z") parameter.dimension[4] = eval.Evaluate(attValue);
   }

   parameter.dimension[0] *= 0.5*lunit;
   parameter.dimension[1] *= 0.5*lunit;
   parameter.dimension[2] *= 0.5*lunit;
   parameter.dimension[3] *= 0.5*lunit;
   parameter.dimension[4] *= 0.5*lunit;
}

void G4GDMLReadParamvol::trap_dimensionsRead(const xercesc::DOMElement* const element,G4GDMLParameterisation::PARAMETER& parameter) {

   G4double lunit = 1.0;
   G4double aunit = 1.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="lunit") lunit = eval.Evaluate(attValue); else
      if (attName=="aunit") aunit = eval.Evaluate(attValue); else
      if (attName=="z") parameter.dimension[0] = eval.Evaluate(attValue); else
      if (attName=="theta") parameter.dimension[1] = eval.Evaluate(attValue); else
      if (attName=="phi") parameter.dimension[2] = eval.Evaluate(attValue); else
      if (attName=="y1") parameter.dimension[3] = eval.Evaluate(attValue); else
      if (attName=="x1") parameter.dimension[4] = eval.Evaluate(attValue); else
      if (attName=="x2") parameter.dimension[5] = eval.Evaluate(attValue); else
      if (attName=="alpha1") parameter.dimension[6] = eval.Evaluate(attValue); else
      if (attName=="y2") parameter.dimension[7] = eval.Evaluate(attValue); else
      if (attName=="x3") parameter.dimension[8] = eval.Evaluate(attValue); else
      if (attName=="x4") parameter.dimension[9] = eval.Evaluate(attValue); else
      if (attName=="alpha2") parameter.dimension[10] = eval.Evaluate(attValue);
   }

   parameter.dimension[0] *= 0.5*lunit;
   parameter.dimension[1] *= aunit;
   parameter.dimension[2] *= aunit;
   parameter.dimension[3] *= 0.5*lunit;
   parameter.dimension[4] *= 0.5*lunit;
   parameter.dimension[5] *= 0.5*lunit;
   parameter.dimension[6] *= aunit;
   parameter.dimension[7] *= 0.5*lunit;
   parameter.dimension[8] *= 0.5*lunit;
   parameter.dimension[9] *= 0.5*lunit;
   parameter.dimension[10] *= aunit;
}

void G4GDMLReadParamvol::tube_dimensionsRead(const xercesc::DOMElement* const element,G4GDMLParameterisation::PARAMETER& parameter) {

   G4double lunit = 1.0;
   G4double aunit = 1.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="lunit") lunit = eval.Evaluate(attValue); else
      if (attName=="aunit") aunit = eval.Evaluate(attValue); else
      if (attName=="rmin") parameter.dimension[0] = eval.Evaluate(attValue); else
      if (attName=="rmax") parameter.dimension[1] = eval.Evaluate(attValue); else
      if (attName=="z") parameter.dimension[2] = eval.Evaluate(attValue); else
      if (attName=="startphi") parameter.dimension[3] = eval.Evaluate(attValue); else
      if (attName=="deltaphi") parameter.dimension[4] = eval.Evaluate(attValue);
   }

   parameter.dimension[0] *= lunit;
   parameter.dimension[1] *= lunit;
   parameter.dimension[2] *= 0.5*lunit;
   parameter.dimension[3] *= aunit;
   parameter.dimension[4] *= aunit;
}

void G4GDMLReadParamvol::cone_dimensionsRead(const xercesc::DOMElement* const element,G4GDMLParameterisation::PARAMETER& parameter) {

   G4double lunit = 1.0;
   G4double aunit = 1.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="lunit") lunit = eval.Evaluate(attValue); else
      if (attName=="aunit") aunit = eval.Evaluate(attValue); else
      if (attName=="rmin1") parameter.dimension[0] = eval.Evaluate(attValue); else
      if (attName=="rmax1") parameter.dimension[1] = eval.Evaluate(attValue); else
      if (attName=="rmin2") parameter.dimension[2] = eval.Evaluate(attValue); else
      if (attName=="rmax2") parameter.dimension[3] = eval.Evaluate(attValue); else
      if (attName=="z") parameter.dimension[4] = eval.Evaluate(attValue); else
      if (attName=="startphi") parameter.dimension[5] = eval.Evaluate(attValue); else
      if (attName=="deltaphi") parameter.dimension[6] = eval.Evaluate(attValue);
   }

   parameter.dimension[0] *= lunit;
   parameter.dimension[1] *= lunit;
   parameter.dimension[2] *= lunit;
   parameter.dimension[3] *= lunit;
   parameter.dimension[4] *= 0.5*lunit;
   parameter.dimension[5] *= aunit;
   parameter.dimension[6] *= aunit;
}

void G4GDMLReadParamvol::sphere_dimensionsRead(const xercesc::DOMElement* const element,G4GDMLParameterisation::PARAMETER& parameter) {

   G4double lunit = 1.0;
   G4double aunit = 1.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="lunit") lunit = eval.Evaluate(attValue); else
      if (attName=="aunit") aunit = eval.Evaluate(attValue); else
      if (attName=="rmin") parameter.dimension[0] = eval.Evaluate(attValue); else
      if (attName=="rmax") parameter.dimension[1] = eval.Evaluate(attValue); else
      if (attName=="startphi") parameter.dimension[2] = eval.Evaluate(attValue); else
      if (attName=="deltaphi") parameter.dimension[3] = eval.Evaluate(attValue); else
      if (attName=="starttheta") parameter.dimension[4] = eval.Evaluate(attValue); else
      if (attName=="deltatheta") parameter.dimension[5] = eval.Evaluate(attValue);
   }

   parameter.dimension[0] *= lunit;
   parameter.dimension[1] *= lunit;
   parameter.dimension[2] *= aunit;
   parameter.dimension[3] *= aunit;
   parameter.dimension[4] *= aunit;
   parameter.dimension[5] *= aunit;
}

void G4GDMLReadParamvol::orb_dimensionsRead(const xercesc::DOMElement* const element,G4GDMLParameterisation::PARAMETER& parameter) {

   G4double lunit = 1.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="lunit") lunit = eval.Evaluate(attValue); else
      if (attName=="r") parameter.dimension[0] = eval.Evaluate(attValue);
   }

   parameter.dimension[0] *= lunit;
}

void G4GDMLReadParamvol::torus_dimensionsRead(const xercesc::DOMElement* const element,G4GDMLParameterisation::PARAMETER& parameter) {

   G4double lunit = 1.0;
   G4double aunit = 1.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="lunit") lunit = eval.Evaluate(attValue); else
      if (attName=="aunit") aunit = eval.Evaluate(attValue); else
      if (attName=="rmin") parameter.dimension[0] = eval.Evaluate(attValue); else
      if (attName=="rmax") parameter.dimension[1] = eval.Evaluate(attValue); else
      if (attName=="rtor") parameter.dimension[2] = eval.Evaluate(attValue); else
      if (attName=="startphi") parameter.dimension[3] = eval.Evaluate(attValue); else
      if (attName=="deltaphi") parameter.dimension[4] = eval.Evaluate(attValue);
   }

   parameter.dimension[0] *= lunit;
   parameter.dimension[1] *= lunit;
   parameter.dimension[2] *= lunit;
   parameter.dimension[3] *= aunit;
   parameter.dimension[4] *= aunit;
}

void G4GDMLReadParamvol::para_dimensionsRead(const xercesc::DOMElement* const element,G4GDMLParameterisation::PARAMETER& parameter) {

   G4double lunit = 1.0;
   G4double aunit = 1.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="lunit") lunit = eval.Evaluate(attValue); else
      if (attName=="aunit") aunit = eval.Evaluate(attValue); else
      if (attName=="x") parameter.dimension[0] = eval.Evaluate(attValue); else
      if (attName=="y") parameter.dimension[1] = eval.Evaluate(attValue); else
      if (attName=="z") parameter.dimension[2] = eval.Evaluate(attValue); else
      if (attName=="alpha") parameter.dimension[3] = eval.Evaluate(attValue); else
      if (attName=="theta") parameter.dimension[4] = eval.Evaluate(attValue); else
      if (attName=="phi") parameter.dimension[5] = eval.Evaluate(attValue);
   }

   parameter.dimension[0] = 0.5*lunit;
   parameter.dimension[1] = 0.5*lunit;
   parameter.dimension[2] = 0.5*lunit;
   parameter.dimension[3] = aunit;
   parameter.dimension[4] = aunit;
   parameter.dimension[5] = aunit;
}

void G4GDMLReadParamvol::hype_dimensionsRead(const xercesc::DOMElement* const element,G4GDMLParameterisation::PARAMETER& parameter) {

   G4double lunit = 1.0;
   G4double aunit = 1.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="lunit") lunit = eval.Evaluate(attValue); else
      if (attName=="aunit") aunit = eval.Evaluate(attValue); else
      if (attName=="rmin") parameter.dimension[0] = eval.Evaluate(attValue); else
      if (attName=="rmax") parameter.dimension[1] = eval.Evaluate(attValue); else
      if (attName=="inst") parameter.dimension[2] = eval.Evaluate(attValue); else
      if (attName=="outst") parameter.dimension[3] = eval.Evaluate(attValue); else
      if (attName=="z") parameter.dimension[4] = eval.Evaluate(attValue);
   }

   parameter.dimension[0] = lunit;
   parameter.dimension[1] = lunit;
   parameter.dimension[2] = aunit;
   parameter.dimension[3] = aunit;
   parameter.dimension[4] = 0.5*lunit;
}

void G4GDMLReadParamvol::parametersRead(const xercesc::DOMElement* const element) {

   G4ThreeVector rotation(0.0,0.0,0.0);
   G4ThreeVector position(0.0,0.0,0.0);

   G4GDMLParameterisation::PARAMETER parameter;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = Transcode(child->getTagName());
      
      if (tag=="rotation") vectorRead(child,rotation); else
      if (tag=="position") vectorRead(child,position); else
      if (tag=="box_dimensions") box_dimensionsRead(child,parameter); else
      if (tag=="trd_dimensions") trd_dimensionsRead(child,parameter); else
      if (tag=="trap_dimensions") trap_dimensionsRead(child,parameter); else
      if (tag=="tube_dimensions") tube_dimensionsRead(child,parameter); else
      if (tag=="cone_dimensions") cone_dimensionsRead(child,parameter); else
      if (tag=="sphere_dimensions") cone_dimensionsRead(child,parameter); else
      if (tag=="orb_dimensions") cone_dimensionsRead(child,parameter); else
      if (tag=="torus_dimensions") cone_dimensionsRead(child,parameter); else
      if (tag=="para_dimensions") cone_dimensionsRead(child,parameter); else
      if (tag=="hype_dimensions") hype_dimensionsRead(child,parameter); else
      G4Exception("G4GDML: ERROR! Unknown tag in parameters: "+tag);
   }

   parameter.pRot = new G4RotationMatrix();
   
   parameter.pRot->rotateX(rotation.x());
   parameter.pRot->rotateY(rotation.y());
   parameter.pRot->rotateZ(rotation.z());

   parameter.position = position;

   parameterisation->addParameter(parameter);
}

void G4GDMLReadParamvol::paramvol_contentRead(const xercesc::DOMElement* const element) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = Transcode(child->getTagName());

       if (tag=="parameters") parametersRead(child); else
       if (tag=="loop") loopRead(child,&G4GDMLRead::paramvol_contentRead);
    }
}

void G4GDMLReadParamvol::paramvolRead(const xercesc::DOMElement* const element,G4LogicalVolume* mother) {

   G4String volumeref;

   parameterisation = new G4GDMLParameterisation();

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = Transcode(child->getTagName());

       if (tag=="volumeref") volumeref = refRead(child);
   }

   paramvol_contentRead(element);

   G4LogicalVolume* logvol = getVolume(GenerateName(volumeref));

   if (parameterisation->getSize()==0) G4Exception("G4GDML: ERROR! No parameters are defined in parameterised volume!");

   new G4PVParameterised("",logvol,mother,kUndefined,parameterisation->getSize(),parameterisation);
}
