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
// $Id: G4GDMLMaterials.cc,v 1.12 2007/11/30 15:20:00 gcosmo Exp $
// GEANT4 tag $ Name:$
//
// class G4GDMLMaterials Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLMaterials.hh"

G4double G4GDMLMaterials::atomRead(const xercesc::DOMElement* const element) {

   G4String value;
   G4String unit("g/mole");

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="value") value = attribute_value; else
      if (attribute_name=="unit") unit = attribute_value;
   }

   return eval.Evaluate(value)*eval.Evaluate(unit);
}

G4double G4GDMLMaterials::DRead(const xercesc::DOMElement* const element) {

   G4String value;
   G4String unit("g/cm3");

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="value") value = attribute_value; else
      if (attribute_name=="unit") unit = attribute_value;
   }

   return eval.Evaluate(value)*eval.Evaluate(unit);
}

void G4GDMLMaterials::elementRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String formula;
   G4String Z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="formula") formula = attribute_value; else
      if (attribute_name=="Z") Z = attribute_value;
   }

   G4double _Z = eval.Evaluate(Z);
   G4double _a = 0;

   G4int nComponents = 0;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="atom") _a = atomRead(child); else
      if (tag=="fraction") nComponents++; else
      G4Exception("GDML: Unknown tag in element: "+tag);
   }

   if (nComponents>0) mixtureRead(element,new G4Element(GenerateName(name),formula,nComponents));
   else new G4Element(GenerateName(name),formula,_Z,_a);
}

G4double G4GDMLMaterials::fractionRead(const xercesc::DOMElement* const element,G4String& ref) {

   G4String n;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="n") n = attribute_value; else
      if (attribute_name=="ref") ref = attribute_value;
   }

   return eval.Evaluate(n);
}

void G4GDMLMaterials::isotopeRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String Z;
   G4String N;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="Z") Z = attribute_value; else
      if (attribute_name=="N") N = attribute_value;
   }

   G4double __Z = eval.Evaluate(Z);
   G4double __N = eval.Evaluate(N);
   G4double __a = 0;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="atom") __a = atomRead(child);
   }

   new G4Isotope(GenerateName(name),(G4int)__Z,(G4int)__N,__a);
}

void G4GDMLMaterials::materialRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String Z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="Z") Z = attribute_value;
   }
  
   G4double _Z = eval.Evaluate(Z);
   G4double _D = 0.0;
   G4double _a = 0.0;

   G4int nComponents = 0;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="D") _D = DRead(child);    else
      if (tag=="atom") _a = atomRead(child); else
      if (tag=="fraction") nComponents++; else
      G4Exception("GDML: Unknown tag in material: "+tag);
   }

   if (nComponents>0) mixtureRead(element,new G4Material(GenerateName(name),_D,nComponents));
   else new G4Material(GenerateName(name),_Z,_a,_D);
}

void G4GDMLMaterials::mixtureRead(const xercesc::DOMElement *const element,G4Element *ele) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="fraction") {

         G4String ref;
	 G4double _n = fractionRead(child,ref);

         ele->AddIsotope(getIsotope(GenerateName(ref)),_n);
      } else G4Exception("GDML: Unknown tag in mixture element: "+tag);
   }
}

void G4GDMLMaterials::mixtureRead(const xercesc::DOMElement *const element,G4Material *material) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="D") { /*already processed*/ } else
      if (tag=="fraction") {

         G4String ref;
	 G4double _n = fractionRead(child,ref);
	 
         G4Material *materialPtr = G4Material::GetMaterial(GenerateName(ref),false);
         G4Element *elementPtr = G4Element::GetElement(GenerateName(ref),false);

         if (materialPtr != 0) material->AddMaterial(materialPtr,_n); else
	 if (elementPtr != 0) material->AddElement(elementPtr,_n);

         if ((materialPtr == 0) && (elementPtr == 0)) G4Exception("GDML: Referenced material/element '"+GenerateName(ref)+"' was not found!");   
      } else G4Exception("GDML: Unknown tag in mixture material: "+tag);
   }
}

void G4GDMLMaterials::materialsRead(const xercesc::DOMElement* const element) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="element") elementRead(child); else 
      if (tag=="isotope") isotopeRead(child); else 
      if (tag=="material") materialRead(child); else 
      G4Exception("GDML: Unknown tag in materials: "+tag);
   }
}

G4Isotope* G4GDMLMaterials::getIsotope(const G4String& ref) const {

   G4Isotope* isotopePtr = G4Isotope::GetIsotope(ref,false);

   if (!isotopePtr) G4Exception("GDML: Referenced isotope '"+ref+"' was not found!"); 

   return isotopePtr;
}

G4Material* G4GDMLMaterials::getMaterial(const G4String& ref) const {

   G4Material *materialPtr = G4Material::GetMaterial(ref,false);

   if (!materialPtr) G4Exception("GDML: Referenced material '"+ref+"' was not found!"); 

   return materialPtr;
}
