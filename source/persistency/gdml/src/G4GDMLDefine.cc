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
//
// $Id: G4GDMLDefine.cc,v 1.14 2007/11/30 13:27:24 ztorzsok Exp $
// GEANT4 tag $ Name:$
//
// class G4GDMLDefine Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLDefine.hh"

void G4GDMLDefine::constantRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String value;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name  = attribute_value; else
      if (attribute_name=="value") value = attribute_value;
   }

   eval.defineConstant(name,eval.Evaluate(value));
}

void G4GDMLDefine::positionRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String unit("1");
   G4String x;
   G4String y;
   G4String z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="unit") unit = attribute_value; else
      if (attribute_name=="x") x = attribute_value; else
      if (attribute_name=="y") y = attribute_value; else
      if (attribute_name=="z") z = attribute_value;
   }

   G4double _unit = eval.Evaluate(unit);

   G4double _x = eval.Evaluate(x)*_unit;
   G4double _y = eval.Evaluate(y)*_unit;
   G4double _z = eval.Evaluate(z)*_unit;

   positionMap[GenerateName(name)] = new G4ThreeVector(_x,_y,_z);
}

void G4GDMLDefine::rotationRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String unit("1");
   G4String x;
   G4String y;
   G4String z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="unit") unit = attribute_value; else
      if (attribute_name=="x") x = attribute_value; else
      if (attribute_name=="y") y = attribute_value; else
      if (attribute_name=="z") z = attribute_value;
   }

   G4double _unit = eval.Evaluate(unit);

   G4double _x = eval.Evaluate(x)*_unit;
   G4double _y = eval.Evaluate(y)*_unit;
   G4double _z = eval.Evaluate(z)*_unit;

   rotationMap[GenerateName(name)] = new G4ThreeVector(_x,_y,_z);
}

void G4GDMLDefine::scaleRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String x;
   G4String y;
   G4String z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="x") x = attribute_value; else
      if (attribute_name=="y") y = attribute_value; else
      if (attribute_name=="z") z = attribute_value;
   }

   G4double _x = eval.Evaluate(x);
   G4double _y = eval.Evaluate(y);
   G4double _z = eval.Evaluate(z);

   scaleMap[GenerateName(name)] = new G4ThreeVector(_x,_y,_z);
}

void G4GDMLDefine::variableRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String value;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* node = attributes->item(attribute_index);

      if (node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="value") value = attribute_value;
   }

   eval.defineVariable(name,eval.Evaluate(value));
}

void G4GDMLDefine::defineRead(const xercesc::DOMElement* const element) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="constant") constantRead(child); else
      if (tag=="position") positionRead(child); else
      if (tag=="rotation") rotationRead(child); else
      if (tag=="scale") scaleRead(child); else
      if (tag=="variable") variableRead(child); else
      G4Exception("GDML: Unknown tag in define: "+tag);
   }
}

G4ThreeVector* G4GDMLDefine::getPosition(const G4String& ref) {

   if (positionMap.find(ref) != positionMap.end()) return positionMap[ref];

   G4Exception("GDML: Referenced position '"+ref+"' was not found!");

   return 0;
}

G4ThreeVector* G4GDMLDefine::getRotation(const G4String& ref) {

   if (rotationMap.find(ref) != rotationMap.end()) return rotationMap[ref];

   G4Exception("GDML: Referenced rotation '"+ref+"' was not found!");

   return 0;
}

G4ThreeVector* G4GDMLDefine::getScale(const G4String& ref) {

   if (scaleMap.find(ref) != scaleMap.end()) return scaleMap[ref];

   G4Exception("GDML: Referenced scale '"+ref+"' was not found!");

   return 0;
}
