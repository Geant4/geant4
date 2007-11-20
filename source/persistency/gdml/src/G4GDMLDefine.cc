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
// $Id: G4GDMLDefine.cc,v 1.8 2007-11-20 09:37:11 gcosmo Exp $
// GEANT4 tag $ Name:$
//
// class G4GDMLDefine Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLDefine.hh"

G4GDMLDefine::G4GDMLDefine() {
}

G4GDMLDefine::~G4GDMLDefine() {
}

G4bool G4GDMLDefine::constantRead(const xercesc::DOMElement* const element) {

   G4String name,value;

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

G4bool G4GDMLDefine::positionRead(const xercesc::DOMElement* const element) {

   G4String name,unit,x,y,z;

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

   G4double _x,_y,_z;

   if (!evaluator->Evaluate(_x,x,unit)) return false;
   if (!evaluator->Evaluate(_y,y,unit)) return false;
   if (!evaluator->Evaluate(_z,z,unit)) return false;

   positionMap[prename+name] = new G4ThreeVector(_x,_y,_z);

   return true;
}

G4bool G4GDMLDefine::rotationRead(const xercesc::DOMElement* const element) {

   G4String name,unit,x,y,z;

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
   
   G4double _x,_y,_z;

   if (!evaluator->Evaluate(_x,x,unit)) return false;
   if (!evaluator->Evaluate(_y,y,unit)) return false;
   if (!evaluator->Evaluate(_z,z,unit)) return false;

   rotationMap[prename+name] = new G4ThreeVector(_x,_y,_z);

   return true;
}

G4bool G4GDMLDefine::Read(const xercesc::DOMElement* const element,G4GDMLEvaluator *eval,const G4String& module) {

   evaluator = eval;
   prename = module;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

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

G4ThreeVector* G4GDMLDefine::getPosition(const G4String& ref) {

   if (positionMap.find(ref) != positionMap.end()) return positionMap[ref];

   G4cout << "GDML: Error! Referenced position '" << ref << "' was not found!" << G4endl;

   return 0;
}

G4ThreeVector* G4GDMLDefine::getRotation(const G4String& ref) {

   if (rotationMap.find(ref) != rotationMap.end()) return rotationMap[ref];

   G4cout << "GDML: Error! Referenced rotation '" << ref << "' was not found!" << G4endl;

   return 0;
}
