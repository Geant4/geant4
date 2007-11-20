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
// $Id: G4GDMLMaterials.cc,v 1.7 2007-11-20 09:37:11 gcosmo Exp $
// GEANT4 tag $ Name:$
//
// class G4GDMLMaterials Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLMaterials.hh"

G4bool G4GDMLMaterials::atomRead(const xercesc::DOMElement* const element,G4double& _value) {

   G4String value;
   G4String unit = "g/mole";

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="value") { value = attribute_value; } else
      if (attribute_name=="unit ") { unit  = attribute_value; }
   }

   return evaluator->Evaluate(_value,value,unit);
}

G4bool G4GDMLMaterials::DRead(const xercesc::DOMElement* const element,G4double& _value) {

   G4String value;
   G4String unit = "g/cm3";

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="value") { value = attribute_value; } else
      if (attribute_name=="unit ") { unit  = attribute_value; }
   }

   return evaluator->Evaluate(_value,value,unit);
}

G4bool G4GDMLMaterials::elementRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String formula;
   G4String Z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name"   ) { name    = attribute_value; } else
      if (attribute_name=="formula") { formula = attribute_value; } else
      if (attribute_name=="Z"      ) { Z       = attribute_value; }
   }

   G4double _Z;
   G4double _a;

   if (!evaluator->Evaluate(_Z,Z)) return false;

   G4int nComponents = 0;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="atom"    ) { if (!atomRead(child,_a)) return false; } else
      if (tag=="fraction") { nComponents++;                         } else
      {
         std::cout << "ERROR! Unsupported tag in element: " << tag << std::endl;
         return false;
      }
   }

   if (nComponents > 0) return mixtureRead(element,new G4Element(prename+name,formula,nComponents));

   new G4Element(prename+name,formula,_Z,_a);

   return true;
}

G4bool G4GDMLMaterials::fractionRead(const xercesc::DOMElement* const element,G4double& _n,G4String& ref) {

   G4String n;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="n"  ) { n   = attribute_value; } else
      if (attribute_name=="ref") { ref = attribute_value; }
   }

   return evaluator->Evaluate(_n,n);
}

G4bool G4GDMLMaterials::isotopeRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String Z;
   G4String N;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") { name = attribute_value; } else
      if (attribute_name=="Z"   ) { Z    = attribute_value; } else
      if (attribute_name=="N"   ) { N    = attribute_value; }
   }

   G4double _ZZ, _NN, _aa;

   if (!evaluator->Evaluate(_ZZ,Z)) return false;
   if (!evaluator->Evaluate(_NN,N)) return false;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="atom") { if (!atomRead(child,_aa)) return false; }
   }

   new G4Isotope(prename+name,(G4int)_ZZ,(G4int)_NN,_aa);

   return true;
}

G4bool G4GDMLMaterials::materialRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String Z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") { name = attribute_value; } else
      if (attribute_name=="Z"   ) { Z    = attribute_value; }
   }
  
   G4double _Z;
   G4double _D;
   G4double _a;

   if (!evaluator->Evaluate(_Z,Z)) return false;

   G4int nComponents = 0;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="D"       ) { if (!DRead   (child,_D)) return false; } else
      if (tag=="atom"    ) { if (!atomRead(child,_a)) return false; } else
      if (tag=="fraction") { nComponents++;                         } else
      {
	 G4cout << "GDML: Error! Unknown tag in material: " << tag << G4endl;
         return false;
      }
   }

   if (nComponents > 0) return mixtureRead(element,new G4Material(prename+name,_D,nComponents));

   new G4Material(prename+name,_Z,_a,_D);

   return true;
}

G4bool G4GDMLMaterials::mixtureRead(const xercesc::DOMElement *const element,G4Element *ele) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="fraction") {

         G4String ref;
	 double _n;
	 
	 ref = prename + ref;
	 
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

G4bool G4GDMLMaterials::mixtureRead(const xercesc::DOMElement *const element,G4Material *material) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="D"       ) { /*already processed*/ } else
      if (tag=="fraction") {

         G4String ref;
	 double _n;
	 
	 if (!fractionRead(child,_n,ref)) return false;

         ref = prename + ref;

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

G4bool G4GDMLMaterials::Read(const xercesc::DOMElement* const element,G4GDMLEvaluator* eval,const G4String& module) {

   evaluator = eval;
   prename = module;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="element" ) { if (!elementRead (child)) return false; } else 
      if (tag=="isotope" ) { if (!isotopeRead (child)) return false; } else 
      if (tag=="material") { if (!materialRead(child)) return false; } else 
      {
	 G4cout << "GDML: Error! Unknown tag in materials: " << tag << G4endl;
         return false;
      }
   }

   return true;
}

G4Material* G4GDMLMaterials::getMaterial(const G4String& ref) const {

   G4Material *materialPtr = G4Material::GetMaterial(ref,false);

   if (!materialPtr) G4cout << "GDML: Error! Referenced material '" << ref << "' was not found!" << G4endl;   

   return materialPtr;
}
