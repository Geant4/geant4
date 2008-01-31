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

#include "G4GDMLReadDefine.hh"

void G4GDMLReadDefine::constantRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4double value = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="name") name = attValue; else
      if (attName=="value") value = eval.Evaluate(attValue);
   }

   eval.defineConstant(name,value);
}

void G4GDMLReadDefine::matrixRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String rows;
   G4String cols;
   G4String values;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name  = attribute_value; else
      if (attribute_name=="rows") rows = attribute_value; else
      if (attribute_name=="cols") cols = attribute_value;
   }

   G4int _rows = eval.EvaluateInteger(rows);
   G4int _cols = eval.EvaluateInteger(cols);

   _rows = 0;
   _cols = 0;

   values = xercesc::XMLString::transcode(element->getTextContent());

   G4cout << "Matrix values: " << values << G4endl;
}

void G4GDMLReadDefine::positionRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4double unit = 1.0;
   G4double x = 0.0;
   G4double y = 0.0;
   G4double z = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="name") name = GenerateName(attValue); else
      if (attName=="unit") unit = eval.Evaluate(attValue); else
      if (attName=="x") x = eval.Evaluate(attValue); else
      if (attName=="y") y = eval.Evaluate(attValue); else
      if (attName=="z") z = eval.Evaluate(attValue);
   }

   positionMap[name] = new G4ThreeVector(x*unit,y*unit,z*unit);
}

void G4GDMLReadDefine::rotationRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4double unit = 1.0;
   G4double x = 0.0;
   G4double y = 0.0;
   G4double z = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="name") name = GenerateName(attValue); else
      if (attName=="unit") unit = eval.Evaluate(attValue); else
      if (attName=="x") x = eval.Evaluate(attValue); else
      if (attName=="y") y = eval.Evaluate(attValue); else
      if (attName=="z") z = eval.Evaluate(attValue);
   }

   rotationMap[name] = new G4ThreeVector(x*unit,y*unit,z*unit);
}

void G4GDMLReadDefine::scaleRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4double x = 1.0;
   G4double y = 1.0;
   G4double z = 1.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="name") name = GenerateName(attValue); else
      if (attName=="x") x = eval.Evaluate(attValue); else
      if (attName=="y") y = eval.Evaluate(attValue); else
      if (attName=="z") z = eval.Evaluate(attValue);
   }

   scaleMap[name] = new G4ThreeVector(x,y,z);
}

void G4GDMLReadDefine::variableRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4double value = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="name") name = attValue; else
      if (attName=="value") value = eval.Evaluate(attValue);
   }

   eval.defineVariable(name,value);
}

void G4GDMLReadDefine::quantityRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4double unit = 1.0;
   G4double value = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="name") name = attValue; else
      if (attName=="value") value = eval.Evaluate(attValue); else
      if (attName=="unit") unit = eval.Evaluate(attValue);
   }

   quantityMap[name] = value*unit;
}

void G4GDMLReadDefine::defineRead(const xercesc::DOMElement* const element) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="constant") constantRead(child); else
      if (tag=="matrix") matrixRead(child); else
      if (tag=="position") positionRead(child); else
      if (tag=="rotation") rotationRead(child); else
      if (tag=="scale") scaleRead(child); else
      if (tag=="variable") variableRead(child); else
      if (tag=="quantity") quantityRead(child); else
      G4Exception("GDML Reader: ERROR! Unknown tag in define: "+tag);
   }
}

G4ThreeVector* G4GDMLReadDefine::getPosition(const G4String& ref) {

   if (positionMap.find(ref) == positionMap.end()) G4Exception("GDML Reader: ERROR! Referenced position '"+ref+"' was not found!");

   return positionMap[ref];
}

G4ThreeVector* G4GDMLReadDefine::getRotation(const G4String& ref) {

   if (rotationMap.find(ref) == rotationMap.end()) G4Exception("GDML Reader: ERROR! Referenced rotation '"+ref+"' was not found!");

   return rotationMap[ref];
}

G4ThreeVector* G4GDMLReadDefine::getScale(const G4String& ref) {

   if (scaleMap.find(ref) == scaleMap.end()) G4Exception("GDML Reader: ERROR! Referenced scale '"+ref+"' was not found!");

   return scaleMap[ref];
}

G4double G4GDMLReadDefine::getQuantity(const G4String& ref) {

   if (quantityMap.find(ref) == quantityMap.end()) G4Exception("GDML Reader: ERROR! Referenced quantity '"+ref+"' was not found!");

   return quantityMap[ref];
}

