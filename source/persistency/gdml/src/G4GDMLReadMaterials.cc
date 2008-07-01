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
// $Id: G4GDMLReadMaterials.cc,v 1.11 2008-07-01 08:12:32 gcosmo Exp $
// GEANT4 tag $ Name:$
//
// class G4GDMLReadMaterials Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLReadMaterials.hh"

G4double G4GDMLReadMaterials::atomRead(const xercesc::DOMElement* const atomElement) {

   G4double value = 0.0;
   G4double unit = g/mole;

   const xercesc::DOMNamedNodeMap* const attributes = atomElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="value") value = eval.Evaluate(attValue); else
      if (attName=="unit") unit = eval.Evaluate(attValue);
   }

   return value*unit;
}

G4int G4GDMLReadMaterials::compositeRead(const xercesc::DOMElement* const compositeElement,G4String& ref) {

   G4int n = 0;

   const xercesc::DOMNamedNodeMap* const attributes = compositeElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="n") n = eval.EvaluateInteger(attValue); else
      if (attName=="ref") ref = attValue;
   }

   return n;
}

G4double G4GDMLReadMaterials::DRead(const xercesc::DOMElement* const DElement) {

   G4double value = 0.0;
   G4double unit = g/cm3;

   const xercesc::DOMNamedNodeMap* const attributes = DElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="value") value = eval.Evaluate(attValue); else
      if (attName=="unit") unit = eval.Evaluate(attValue);
   }

   return value*unit;
}

G4double G4GDMLReadMaterials::PRead(const xercesc::DOMElement* const PElement) {

   G4double value = STP_Pressure;
   G4double unit = pascal;

   const xercesc::DOMNamedNodeMap* const attributes = PElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="value") value = eval.Evaluate(attValue); else
      if (attName=="unit") unit = eval.Evaluate(attValue);
   }

   return value*unit;
}

G4double G4GDMLReadMaterials::TRead(const xercesc::DOMElement* const TElement) {

   G4double value = STP_Temperature;
   G4double unit = kelvin;

   const xercesc::DOMNamedNodeMap* const attributes = TElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="value") value = eval.Evaluate(attValue); else
      if (attName=="unit") unit = eval.Evaluate(attValue);
   }

   return value*unit;
}

void G4GDMLReadMaterials::elementRead(const xercesc::DOMElement* const elementElement) {

   G4String name;
   G4String formula;
   G4double a = 0.0;
   G4double Z = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes = elementElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name") name = GenerateName(attValue); else
      if (attName=="formula") formula = attValue; else
      if (attName=="Z") Z = eval.Evaluate(attValue);
   }

   G4int nComponents = 0;

   for (xercesc::DOMNode* iter = elementElement->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = Transcode(child->getTagName());

      if (tag=="atom") a = atomRead(child); else
      if (tag=="fraction") nComponents++;
   }

   if (nComponents>0) mixtureRead(elementElement,new G4Element(name,formula,nComponents));
   else new G4Element(name,formula,Z,a);
}

G4double G4GDMLReadMaterials::fractionRead(const xercesc::DOMElement* const fractionElement,G4String& ref) {

   G4double n = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes = fractionElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="n") n = eval.Evaluate(attValue); else
      if (attName=="ref") ref = attValue;
   }

   return n;
}

void G4GDMLReadMaterials::isotopeRead(const xercesc::DOMElement* const isotopeElement) {

   G4String name;
   G4int Z = 0;
   G4int N = 0;
   G4double a = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes = isotopeElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name") name = GenerateName(attValue); else
      if (attName=="Z") Z = eval.EvaluateInteger(attValue); else
      if (attName=="N") N = eval.EvaluateInteger(attValue);
   }

   for (xercesc::DOMNode* iter = isotopeElement->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = Transcode(child->getTagName());

      if (tag=="atom") a = atomRead(child);
   }

   new G4Isotope(name,Z,N,a);
}

void G4GDMLReadMaterials::materialRead(const xercesc::DOMElement* const materialElement) {

   G4String name;
   G4double Z = 0.0;
   G4double a = 0.0;
   G4double D = 0.0;
   G4State state = kStateUndefined;
   G4double T = STP_Temperature;
   G4double P = STP_Pressure;

   const xercesc::DOMNamedNodeMap* const attributes = materialElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name") name = GenerateName(attValue); else
      if (attName=="Z") Z = eval.Evaluate(attValue); else
      if (attName=="state") {
      
         if (attValue=="solid") state = kStateSolid; else
         if (attValue=="liquid") state = kStateLiquid; else
         if (attValue=="gas") state = kStateGas;
      }
   }

   size_t nComponents = 0;

   for (xercesc::DOMNode* iter = materialElement->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = Transcode(child->getTagName());

      if (tag=="atom") a = atomRead(child); else
      if (tag=="Dref") D = getQuantity(GenerateName(refRead(child))); else
      if (tag=="Pref") P = getQuantity(GenerateName(refRead(child))); else
      if (tag=="Tref") T = getQuantity(GenerateName(refRead(child))); else
      if (tag=="D") D = DRead(child); else
      if (tag=="P") P = PRead(child); else
      if (tag=="T") T = TRead(child); else
      if (tag=="fraction" || tag=="composite") nComponents++;
   }

   G4Material* material =  0;

   if (nComponents==0) material = new G4Material(name,Z,a,D,state,T,P);
   else mixtureRead(materialElement,material = new G4Material(name,D,nComponents,state,T,P));

   for (xercesc::DOMNode* iter = materialElement->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = Transcode(child->getTagName());

      if (tag=="property") propertyRead(child,material);
   }
}

void G4GDMLReadMaterials::mixtureRead(const xercesc::DOMElement *const mixtureElement,G4Element *element) {

   for (xercesc::DOMNode* iter = mixtureElement->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = Transcode(child->getTagName());

      if (tag=="fraction") {

         G4String ref;
	 G4double n = fractionRead(child,ref);
         element->AddIsotope(getIsotope(GenerateName(ref)),n);
      }
   }
}

void G4GDMLReadMaterials::mixtureRead(const xercesc::DOMElement *const mixtureElement,G4Material *material) {

   for (xercesc::DOMNode* iter = mixtureElement->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = Transcode(child->getTagName());

      if (tag=="fraction") {

         G4String ref;
	 G4double n = fractionRead(child,ref);
	 
         G4Material *materialPtr = getMaterial(GenerateName(ref),false);
         G4Element *elementPtr = getElement(GenerateName(ref),false);

         if (materialPtr != 0) material->AddMaterial(materialPtr,n); else
	 if (elementPtr != 0) material->AddElement(elementPtr,n);

         if ((materialPtr == 0) && (elementPtr == 0)) G4Exception("G4GDML: ERROR! Referenced material/element '"+GenerateName(ref)+"' was not found!");   
      } 
      else if (tag=="composite") {
      
         G4String ref;
	 G4int n = compositeRead(child,ref);

         G4Element *elementPtr = getElement(GenerateName(ref));

         material->AddElement(elementPtr,n);
      }
   }
}

void G4GDMLReadMaterials::propertyRead(const xercesc::DOMElement* const propertyElement,G4Material* material) {

   G4String name;
   G4String ref;
   G4GDMLMatrix matrix;

   const xercesc::DOMNamedNodeMap* const attributes = propertyElement->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   
      const G4String attName = Transcode(attribute->getName());
      const G4String attValue = Transcode(attribute->getValue());

      if (attName=="name") name = GenerateName(attValue); else
      if (attName=="ref") matrix = getMatrix(ref=attValue);
   }

   if (matrix.getCols() != 2) G4Exception("G4GDML: ERROR! Referenced matrix '"+ref+"' should have two columns as a property table for material: "+material->GetName());
   if (matrix.getRows() == 0) return;

   G4MaterialPropertiesTable* matprop = material->GetMaterialPropertiesTable();
   if (!matprop) material->SetMaterialPropertiesTable(matprop = new G4MaterialPropertiesTable());

   G4MaterialPropertyVector* propvect = new G4MaterialPropertyVector(0,0,0);
   for (size_t i=0;i<matrix.getRows();i++) propvect->AddElement(matrix.get(i,0),matrix.get(i,1));
   matprop->AddProperty(name,propvect);
}

void G4GDMLReadMaterials::materialsRead(const xercesc::DOMElement* const materialsElement) {

   G4cout << "G4GDML: Reading materials..." << G4endl;

   for (xercesc::DOMNode* iter = materialsElement->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
      const G4String tag = Transcode(child->getTagName());

      if (tag=="element") elementRead(child); else 
      if (tag=="isotope") isotopeRead(child); else 
      if (tag=="material") materialRead(child); else 
      G4Exception("G4GDML: ERROR! Unknown tag in materials: "+tag);
   }
}

G4Element* G4GDMLReadMaterials::getElement(const G4String& ref,bool verbose) const {

   G4Element* elementPtr = G4Element::GetElement(ref,false);

   if (!elementPtr) elementPtr = G4NistManager::Instance()->FindOrBuildElement(ref);

   if (verbose && !elementPtr) G4Exception("G4GDML: ERROR! Referenced element '"+ref+"' was not found!"); 

   return elementPtr;
}

G4Isotope* G4GDMLReadMaterials::getIsotope(const G4String& ref,bool verbose) const {

   G4Isotope* isotopePtr = G4Isotope::GetIsotope(ref,false);

   if (verbose && !isotopePtr) G4Exception("G4GDML: ERROR! Referenced isotope '"+ref+"' was not found!"); 

   return isotopePtr;
}

G4Material* G4GDMLReadMaterials::getMaterial(const G4String& ref,bool verbose) const {

   G4Material *materialPtr = G4Material::GetMaterial(ref,false);

   if (!materialPtr) materialPtr = G4NistManager::Instance()->FindOrBuildMaterial(ref);

   if (verbose && !materialPtr) G4Exception("G4GDML: ERROR! Referenced material '"+ref+"' was not found!"); 

   return materialPtr;
}
