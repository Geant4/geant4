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
// $Id: G4GDMLReadMaterials.cc,v 1.2 2008-01-30 12:38:42 ztorzsok Exp $
// GEANT4 tag $ Name:$
//
// class G4GDMLMaterials Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLReadMaterials.hh"

G4double G4GDMLReadMaterials::atomRead(const xercesc::DOMElement* const element) {

   G4double value = 0.0;
   G4String unit("g/mole");

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="value") value = eval.Evaluate(attValue); else
      if (attName=="unit") unit = attValue;
   }

   return value*eval.Evaluate(unit);
}

G4int G4GDMLReadMaterials::compositeRead(const xercesc::DOMElement* const element,G4String& ref) {

   G4int n = 0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="n") n = eval.EvaluateInteger(attValue); else
      if (attName=="ref") ref = attValue;
   }

   return n;
}

G4double G4GDMLReadMaterials::DRead(const xercesc::DOMElement* const element) {

   G4double value = 0.0;
   G4String unit("g/cm3");

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="value") value = eval.Evaluate(attValue); else
      if (attName=="unit") unit = attValue;
   }

   return value*eval.Evaluate(unit);
}

G4double G4GDMLReadMaterials::PRead(const xercesc::DOMElement* const element) {

   G4double value = STP_Pressure;
   G4String unit("pascal");

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="value") value = eval.Evaluate(attValue); else
      if (attName=="unit") unit = attValue;
   }

   return value*eval.Evaluate(unit);
}

G4double G4GDMLReadMaterials::TRead(const xercesc::DOMElement* const element) {

   G4double value = STP_Temperature;
   G4String unit("K");

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="value") value = eval.Evaluate(attValue); else
      if (attName=="unit") unit = attValue;
   }

   return value*eval.Evaluate(unit);
}

void G4GDMLReadMaterials::elementRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String formula;
   G4double a = 0.0;
   G4double Z = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="name") name = attValue; else
      if (attName=="formula") formula = attValue; else
      if (attName=="Z") Z = eval.Evaluate(attValue);
   }

   G4int nComponents = 0;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="atom") a = atomRead(child); else
      if (tag=="fraction") nComponents++; else
      G4Exception("GDML: Unknown tag in element: "+tag);
   }

   if (nComponents>0) mixtureRead(element,new G4Element(GenerateName(name),formula,nComponents));
   else new G4Element(GenerateName(name),formula,Z,a);
}

G4double G4GDMLReadMaterials::fractionRead(const xercesc::DOMElement* const element,G4String& ref) {

   G4double n = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="n") n = eval.Evaluate(attValue); else
      if (attName=="ref") ref = attValue;
   }

   return n;
}

void G4GDMLReadMaterials::isotopeRead(const xercesc::DOMElement* const element) {

   G4String name;

   G4int Z = 0;
   G4int N = 0;
   G4double a = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="name") name = attValue; else
      if (attName=="Z") Z = eval.EvaluateInteger(attValue); else
      if (attName=="N") N = eval.EvaluateInteger(attValue);
   }

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="atom") a = atomRead(child);
   }

   new G4Isotope(GenerateName(name),Z,N,a);
}

void G4GDMLReadMaterials::materialRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4double Z = 0.0;
   G4double a = 0.0;
   G4double D = 0.0;
   G4State state = kStateUndefined;
   G4double T = STP_Temperature;
   G4double P = STP_Pressure;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="name") name = attValue; else
      if (attName=="Z") Z = eval.Evaluate(attValue); else
      if (attName=="state") {
      
         if (attValue=="solid") state = kStateSolid; else
         if (attValue=="liquid") state = kStateLiquid; else
         if (attValue=="gas") state = kStateGas;
      }
   }

   size_t nComponents = 0;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="D") D = DRead(child); else
      if (tag=="P") P = PRead(child); else
      if (tag=="T") T = TRead(child); else
      if (tag=="atom") a = atomRead(child); else
      if (tag=="fraction") nComponents++; else
      if (tag=="composite") nComponents++;
   }

   if (nComponents==0) new G4Material(GenerateName(name),Z,a,D,state,T,P);
   else mixtureRead(element,new G4Material(GenerateName(name),D,nComponents,state,T,P));
}

void G4GDMLReadMaterials::mixtureRead(const xercesc::DOMElement *const element,G4Element *ele) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="fraction") {

         G4String ref;
	 G4double n = fractionRead(child,ref);

         ele->AddIsotope(getIsotope(GenerateName(ref)),n);
      }
   }
}

void G4GDMLReadMaterials::mixtureRead(const xercesc::DOMElement *const element,G4Material *material) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="D") { /*already processed*/ } else
      if (tag=="P") { /*already processed*/ } else
      if (tag=="T") { /*already processed*/ } else
      if (tag=="fraction") {

         G4String ref;
	 G4double n = fractionRead(child,ref);
	 
         G4Material *materialPtr = getMaterial(GenerateName(ref),false);
         G4Element *elementPtr = getElement(GenerateName(ref),false);

         if (materialPtr != 0) material->AddMaterial(materialPtr,n); else
	 if (elementPtr != 0) material->AddElement(elementPtr,n);

         if ((materialPtr == 0) && (elementPtr == 0)) G4Exception("GDML: Referenced material/element '"+GenerateName(ref)+"' was not found!");   
      } 
      else if (tag=="composite") {
      
         G4String ref;
	 G4int n = compositeRead(child,ref);

         G4Element *elementPtr = getElement(GenerateName(ref));

         material->AddElement(elementPtr,n);
      }
   }
}

void G4GDMLReadMaterials::opticalsurfaceRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String smodel;
   G4String sfinish;
   G4String stype;
   G4double value = 0.0;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attName = xercesc::XMLString::transcode(attribute->getName());
      const G4String attValue = xercesc::XMLString::transcode(attribute->getValue());

      if (attName=="name") name = GenerateName(attValue); else
      if (attName=="model") smodel = attValue; else
      if (attName=="finish") sfinish = attValue; else
      if (attName=="type") stype = attValue; else
      if (attName=="value") value = eval.Evaluate(attValue);
   }

   G4OpticalSurfaceModel model; 
   G4OpticalSurfaceFinish finish;
   G4SurfaceType type;   
   
   if (smodel="unified") model = unified; else 
   model = glisur;

   if (sfinish=="polishedfrontpainted") finish = polishedfrontpainted; else
   if (sfinish=="polishedbackpainted") finish = polishedbackpainted; else
   if (sfinish=="groundfrontpainted") finish = groundfrontpainted; else
   if (sfinish=="groundbackpainted") finish = groundbackpainted; else
   if (sfinish=="ground") finish = ground; else
   finish = polished;

   if (stype=="dielectric_metal") type = dielectric_metal; else
   if (stype=="x_ray") type = x_ray; else
   if (stype=="firsov") type = firsov; else   
   type = dielectric_dielectric;

   new G4OpticalSurface(name,model,finish,type,value);
}

void G4GDMLReadMaterials::materialsRead(const xercesc::DOMElement* const element) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="element") elementRead(child); else 
      if (tag=="isotope") isotopeRead(child); else 
      if (tag=="material") materialRead(child); else 
      if (tag=="opticalsurface") opticalsurfaceRead(child); else 
      G4Exception("GDML: Unknown tag in materials: "+tag);
   }
}

G4Element* G4GDMLReadMaterials::getElement(const G4String& ref,bool verbose) const {

   G4Element* elementPtr = G4Element::GetElement(ref,false);

   if (!elementPtr) elementPtr = G4NistManager::Instance()->FindOrBuildElement(ref);

   if (verbose && !elementPtr) G4Exception("GDML: Referenced element '"+ref+"' was not found!"); 

   return elementPtr;
}

G4Isotope* G4GDMLReadMaterials::getIsotope(const G4String& ref,bool verbose) const {

   G4Isotope* isotopePtr = G4Isotope::GetIsotope(ref,false);

   if (verbose && !isotopePtr) G4Exception("GDML: Referenced isotope '"+ref+"' was not found!"); 

   return isotopePtr;
}

G4Material* G4GDMLReadMaterials::getMaterial(const G4String& ref,bool verbose) const {

   G4Material *materialPtr = G4Material::GetMaterial(ref,false);

   if (!materialPtr) materialPtr = G4NistManager::Instance()->FindOrBuildMaterial(ref);

   if (verbose && !materialPtr) G4Exception("GDML: Referenced material '"+ref+"' was not found!"); 

   return materialPtr;
}
