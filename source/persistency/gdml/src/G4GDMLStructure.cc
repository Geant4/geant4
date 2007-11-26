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
// $Id: G4GDMLStructure.cc,v 1.19 2007-11-26 14:31:32 ztorzsok Exp $
// GEANT4 tag $ Name:$
//
// class G4GDMLStructure Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLStructure.hh"

G4GDMLStructure::G4GDMLStructure() {

   evaluator = new G4GDMLEvaluator();

   parser = 0;
}

EAxis G4GDMLStructure::directionRead(const xercesc::DOMElement* const element) {

   G4String x;
   G4String y;
   G4String z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="x") { x = attribute_value; } else
      if (attribute_name=="y") { y = attribute_value; } else
      if (attribute_name=="z") { z = attribute_value; }
   }

   G4double _x = evaluator->Evaluate(x);
   G4double _y = evaluator->Evaluate(y);
   G4double _z = evaluator->Evaluate(z);

   if (_x == 1.0 && _y == 0.0 && _z == 0.0) return kXAxis; else
   if (_x == 0.0 && _y == 1.0 && _z == 0.0) return kYAxis; else
   if (_x == 0.0 && _y == 0.0 && _z == 1.0) return kZAxis;

   G4Exception("GDML: Only directions along axes are supported!");

   return kZAxis;
}

void G4GDMLStructure::divisionvolRead(const xercesc::DOMElement* const element,G4LogicalVolume* pMother) {

   G4String unit;
   G4String axis;
   G4String width;
   G4String offset;
   G4String number;
   G4String volumeref;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="unit"  ) unit   = attribute_value; else
      if (attribute_name=="axis"  ) axis   = attribute_value; else
      if (attribute_name=="width" ) width  = attribute_value; else
      if (attribute_name=="offset") offset = attribute_value; else
      if (attribute_name=="number") number = attribute_value;
   }

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="volumeref") volumeref = refRead(child);
   }

   G4LogicalVolume* pLogical = getVolume(file+volumeref);

   G4double _unit = evaluator->Evaluate(unit);

   G4double _width  = evaluator->Evaluate(width )*_unit;
   G4double _offset = evaluator->Evaluate(offset)*_unit;
   G4double _number = evaluator->Evaluate(number);

   EAxis    _axis   = kZAxis;

   if (axis=="kXAxis") _axis = kXAxis; else
   if (axis=="kYAxis") _axis = kYAxis; else
   if (axis=="kZAxis") _axis = kZAxis; else
   if (axis=="kRho"  ) _axis = kRho;   else
   if (axis=="kPhi"  ) _axis = kPhi;
   
   new G4PVDivision("",pLogical,pMother,_axis,(G4int)_number,_width,_offset);
}

G4LogicalVolume* G4GDMLStructure::fileRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String volname;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name"   ) name    = attribute_value; else
      if (attribute_name=="volname") volname = attribute_value;
   }

   G4GDMLStructure structure;
   
   structure.gdmlRead(name,parser);

   return structure.getVolume(structure.file + volname);
}

void G4GDMLStructure::loopRead(const xercesc::DOMElement* const element) {

   G4String var;
   G4String from;
   G4String to;
   G4String step;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="var" ) var  = attribute_value; else
      if (attribute_name=="from") from = attribute_value; else
      if (attribute_name=="to"  ) to   = attribute_value; else
      if (attribute_name=="step") step = attribute_value;
   }

   G4double _var  = evaluator->Evaluate(var );
   G4double _from = evaluator->Evaluate(from);
   G4double _to   = evaluator->Evaluate(to  );
   G4double _step = evaluator->Evaluate(step);
   
   if (!from.empty()) _var = _from;

   while (_var <= _to) {
   
      evaluator->setVariable(var,_var);
      Read(element,file);

      _var += _step;
   }
}

void G4GDMLStructure::paramvolRead(const xercesc::DOMElement* const element,G4LogicalVolume *pMotherLogical) {

   G4String volumeref;
   G4String parameterised_position_size;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="parameterised_position_size") { ; }  else
      if (tag=="volumeref"                  ) volumeref = refRead(child);
   }
}

void G4GDMLStructure::physvolRead(const xercesc::DOMElement* const element,G4LogicalVolume *mother) {

   G4String volumeref;
   G4String positionref;
   G4String rotationref;

   G4LogicalVolume* logvol = 0;
   G4ThreeVector tlate;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="file"       ) logvol = fileRead(child); else
      if (tag=="volumeref"  ) volumeref = refRead(child); else
      if (tag=="position"   ) tlate = positionRead(child); else
      if (tag=="positionref") positionref = refRead(child); else
      if (tag=="rotationref") rotationref = refRead(child);
   }

   if (!volumeref.empty()) logvol = getVolume(solids.nameProcess(volumeref));

   if (!positionref.empty()) tlate = *solids.define.getPosition(file+positionref);

   G4RotationMatrix* pRot=0;

   if (!rotationref.empty()) {

      G4ThreeVector* anglePtr = solids.define.getRotation(file+rotationref);

      pRot = new G4RotationMatrix();

      pRot->rotateX(anglePtr->x());
      pRot->rotateY(anglePtr->y());
      pRot->rotateZ(anglePtr->z());
   }
 
   new G4PVPlacement(pRot,tlate,logvol,"",mother,false,0);
}

G4ThreeVector G4GDMLStructure::positionRead(const xercesc::DOMElement* const element) {

   G4String unit;
   G4String x;
   G4String y;
   G4String z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="unit") { unit = attribute_value; } else
      if (attribute_name=="x"   ) { x    = attribute_value; } else
      if (attribute_name=="y"   ) { y    = attribute_value; } else
      if (attribute_name=="z"   ) { z    = attribute_value; }
   }

   G4double _unit = evaluator->Evaluate(unit);

   G4double _x = evaluator->Evaluate(x)*_unit;
   G4double _y = evaluator->Evaluate(y)*_unit;
   G4double _z = evaluator->Evaluate(z)*_unit;
   
   return G4ThreeVector(_x,_y,_z);
}

G4double G4GDMLStructure::quantityRead(const xercesc::DOMElement* const element) {

   G4String value;
   G4String unit;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="value") { value = attribute_value; }
      if (attribute_name=="unit" ) { unit  = attribute_value; }
   }

   return evaluator->Evaluate(value)*evaluator->Evaluate(unit);
}

G4String G4GDMLStructure::refRead(const xercesc::DOMElement* const element) {

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

void G4GDMLStructure::replicate_along_axisRead(const xercesc::DOMElement* const element,G4double& _width,G4double& _offset,EAxis& _axis) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="width"    ) _width = quantityRead(child); else
      if (tag=="offset"   ) _offset = quantityRead(child); else
      if (tag=="direction") _axis = directionRead(child);
   }
}

void G4GDMLStructure::replicavolRead(const xercesc::DOMElement* const element,G4LogicalVolume* pMother) {

   G4String volumeref;
   G4String numb;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="numb") numb = attribute_value;
   }

   G4double _numb = evaluator->Evaluate(numb);
   G4double _width = 0.0;
   G4double _offset = 0.0;
   EAxis _axis;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="volumeref"           ) volumeref = refRead(child); else
      if (tag=="replicate_along_axis") replicate_along_axisRead(child,_width,_offset,_axis);
   }

   G4LogicalVolume* pLogical = getVolume(file+volumeref);

   new G4PVReplica("",pLogical,pMother,_axis,(G4int)_numb,_width,_offset);
}

void G4GDMLStructure::volumeRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String solidref;
   G4String materialref;

   XMLCh *name_attr = xercesc::XMLString::transcode("name");
   name = xercesc::XMLString::transcode(element->getAttribute(name_attr));
   xercesc::XMLString::release(&name_attr);

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="materialref") materialref = refRead(child); else
      if (tag=="solidref"   ) solidref = refRead(child);
   }

   G4Material* materialPtr = materials.getMaterial(file+materialref); 
   G4VSolid* solidPtr = solids.getSolid(solids.nameProcess(solidref));

   volume_contentRead(element,new G4LogicalVolume(solidPtr,materialPtr,solids.nameProcess(name),0,0,0));
}

void G4GDMLStructure::volume_contentRead(const xercesc::DOMElement* const element,G4LogicalVolume* volumePtr) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
  
      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="loop"       ) volume_loopRead(child,volumePtr); else
      if (tag=="paramvol"   ) paramvolRead(child,volumePtr); else
      if (tag=="physvol"    ) physvolRead(child,volumePtr); else
      if (tag=="replicavol" ) replicavolRead(child,volumePtr); else
      if (tag=="divisionvol") divisionvolRead(child,volumePtr);
   }
}

void G4GDMLStructure::volume_loopRead(const xercesc::DOMElement* const element,G4LogicalVolume* volumePtr) {

   G4String var;
   G4String from;
   G4String to;
   G4String step;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="var" ) var  = attribute_value; else
      if (attribute_name=="from") from = attribute_value; else
      if (attribute_name=="to"  ) to   = attribute_value; else
      if (attribute_name=="step") step = attribute_value;
   }

   G4double _var  = evaluator->Evaluate(var );
   G4double _from = evaluator->Evaluate(from);
   G4double _to   = evaluator->Evaluate(to  );
   G4double _step = evaluator->Evaluate(step);
   
   if (!from.empty()) _var = _from;

   while (_var <= _to) {
   
      evaluator->setVariable(var,_var);
      volume_contentRead(element,volumePtr);

      _var += _step;
   }
}

void G4GDMLStructure::Read(const xercesc::DOMElement* const element,const G4String& file0) {

   file = file0;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="loop"  ) loopRead(child); else
      if (tag=="volume") volumeRead(child); else      
      G4Exception("GDML: Unknown tag in structure: "+tag);
   }
}

void G4GDMLStructure::gdmlRead(const G4String& fileName,xercesc::XercesDOMParser* newParser) {

   file = fileName + "_";

   parser = newParser;

   try {

      parser->parse(fileName.c_str());
   }
   catch (const xercesc::XMLException &e) {
   
      char* message = xercesc::XMLString::transcode(e.getMessage());
      G4cout << "XML: " << message << G4endl;
      xercesc::XMLString::release(&message);
   }
   catch (const xercesc::DOMException &e) {
   
      char* message = xercesc::XMLString::transcode(e.getMessage());
      G4cout << "DOM: " << message << G4endl;
      xercesc::XMLString::release(&message);
   }

   xercesc::DOMDocument* doc = parser->getDocument();

   if (!doc) G4Exception("GDML: Unable to open document: "+fileName);

   xercesc::DOMElement* element = doc->getDocumentElement();

   if (!element) G4Exception("GDML: Empty document!");

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="define"   ) solids.define.Read(child,evaluator,file); else
      if (tag=="materials") materials.Read(child,evaluator,file); else
      if (tag=="solids"   ) solids.Read(child,evaluator,file); else
      if (tag=="setup"    ) setup.Read(child,file); else
      if (tag=="structure") Read(child,file);
   }
}

G4LogicalVolume* G4GDMLStructure::getVolume(const G4String& ref) const {

   G4LogicalVolume *volumePtr = G4LogicalVolumeStore::GetInstance()->GetVolume(ref,false);

   if (!volumePtr) G4Exception("GDML: Referenced volume '"+ref+"' was not found!");

   return volumePtr;
}
