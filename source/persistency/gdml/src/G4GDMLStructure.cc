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
// $Id: G4GDMLStructure.cc,v 1.25 2007/11/30 14:51:20 ztorzsok Exp $
// GEANT4 tag $ Name:$
//
// class G4GDMLStructure Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLStructure.hh"

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

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="x") x = attribute_value; else
      if (attribute_name=="y") y = attribute_value; else
      if (attribute_name=="z") z = attribute_value;
   }

   G4double _x = eval.Evaluate(x);
   G4double _y = eval.Evaluate(y);
   G4double _z = eval.Evaluate(z);

   if (_x == 1.0 && _y == 0.0 && _z == 0.0) return kXAxis; else
   if (_x == 0.0 && _y == 1.0 && _z == 0.0) return kYAxis; else
   if (_x == 0.0 && _y == 0.0 && _z == 1.0) return kZAxis;

   G4Exception("GDML: Only directions along axes are supported!");

   return kZAxis;
}

void G4GDMLStructure::divisionvolRead(const xercesc::DOMElement* const element,G4LogicalVolume* pMother) {

   G4String unit("1");
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

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="unit") unit = attribute_value; else
      if (attribute_name=="axis") axis = attribute_value; else
      if (attribute_name=="width") width  = attribute_value; else
      if (attribute_name=="offset") offset = attribute_value; else
      if (attribute_name=="number") number = attribute_value;
   }

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="volumeref") volumeref = refRead(child);
   }

   G4LogicalVolume* pLogical = getVolume(GenerateName(volumeref));

   G4double _unit = eval.Evaluate(unit);

   G4double _width  = eval.Evaluate(width)*_unit;
   G4double _offset = eval.Evaluate(offset)*_unit;
   G4double _number = eval.Evaluate(number);

   EAxis _axis = kZAxis;

   if (axis=="kXAxis") _axis = kXAxis; else
   if (axis=="kYAxis") _axis = kYAxis; else
   if (axis=="kZAxis") _axis = kZAxis; else
   if (axis=="kRho") _axis = kRho; else
   if (axis=="kPhi") _axis = kPhi;
   
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

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="volname") volname = attribute_value;
   }

   G4GDMLStructure structure; // We create a new structure with a new evaluator
   
   structure.Parse(name);

   return structure.getVolume(structure.GenerateName(volname));
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

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="var") var  = attribute_value; else
      if (attribute_name=="from") from = attribute_value; else
      if (attribute_name=="to") to = attribute_value; else
      if (attribute_name=="step") step = attribute_value;
   }

   eval.checkVariable(var);

   G4int _var = eval.EvaluateInteger(var );
   G4int _from = eval.EvaluateInteger(from);
   G4int _to = eval.EvaluateInteger(to  );
   G4int _step = eval.EvaluateInteger(step);
   
   if (!from.empty()) _var = _from;

   while (_var <= _to) {
   
      eval.setVariable(var,_var);
      structureRead(element);

      _var += _step;
   }
}

void G4GDMLStructure::paramvolRead(const xercesc::DOMElement* const element,G4LogicalVolume *pMotherLogical) {

   pMotherLogical = 0;

   G4String volumeref;
   G4String parameterised_position_size;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="parameterised_position_size") { ; }  else
      if (tag=="volumeref") volumeref = refRead(child);
   }
}

void G4GDMLStructure::physvolRead(const xercesc::DOMElement* const element,G4LogicalVolume *mother) {

   G4String volumeref;
   G4String positionref;
   G4String rotationref;
   G4String scaleref;

   G4LogicalVolume* logvol = 0;

   G4ThreeVector position;
   G4ThreeVector rotation;
   G4ThreeVector scale(1.0,1.0,1.0);

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="file") logvol = fileRead(child); else
      if (tag=="volumeref") volumeref = refRead(child); else
      if (tag=="position") position = positionRead(child); else
      if (tag=="rotation") rotation = rotationRead(child); else
      if (tag=="scale") scale = scaleRead(child); else
      if (tag=="positionref") positionref = refRead(child); else
      if (tag=="rotationref") rotationref = refRead(child); else
      if (tag=="scaleref") scaleref = refRead(child);
   }

   if (!volumeref.empty()) logvol = getVolume(GenerateName(volumeref));

   if (!positionref.empty()) position = *getPosition(GenerateName(positionref));
   if (!rotationref.empty()) rotation = *getRotation(GenerateName(rotationref));
   if (!scaleref.empty()) scale = *getScale(GenerateName(scaleref));

   G4RotationMatrix Rot;

   Rot.rotateX(rotation.x());
   Rot.rotateY(rotation.y());
   Rot.rotateZ(rotation.z());
   
   G4Transform3D transform(Rot.inverse(),position);
   transform = transform*G4Scale3D(scale.x(),scale.y(),scale.z());

   G4ReflectionFactory::Instance()->Place(transform,"",logvol,mother,false,0);
}

G4double G4GDMLStructure::quantityRead(const xercesc::DOMElement* const element) {

   G4String value;
   G4String unit("1");

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="value") value = attribute_value;
      if (attribute_name=="unit") unit = attribute_value;
   }

   return eval.Evaluate(value)*eval.Evaluate(unit);
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

      if (tag=="width") _width = quantityRead(child); else
      if (tag=="offset") _offset = quantityRead(child); else
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

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="numb") numb = attribute_value;
   }

   G4double _numb = eval.Evaluate(numb);
   G4double _width = 0.0;
   G4double _offset = 0.0;
   EAxis _axis;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="volumeref") volumeref = refRead(child); else
      if (tag=="replicate_along_axis") replicate_along_axisRead(child,_width,_offset,_axis);
   }

   G4LogicalVolume* pLogical = getVolume(GenerateName(volumeref));

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
      if (tag=="solidref") solidref = refRead(child);
   }

   G4Material* materialPtr = getMaterial(GenerateName(materialref)); 
   G4VSolid* solidPtr = getSolid(GenerateName(solidref));

   volume_contentRead(element,new G4LogicalVolume(solidPtr,materialPtr,GenerateName(name),0,0,0));
}

void G4GDMLStructure::volume_contentRead(const xercesc::DOMElement* const element,G4LogicalVolume* volumePtr) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
  
      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="loop") volume_loopRead(child,volumePtr); else
      if (tag=="paramvol") paramvolRead(child,volumePtr); else
      if (tag=="physvol") physvolRead(child,volumePtr); else
      if (tag=="replicavol") replicavolRead(child,volumePtr); else
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

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="var") var  = attribute_value; else
      if (attribute_name=="from") from = attribute_value; else
      if (attribute_name=="to") to = attribute_value; else
      if (attribute_name=="step") step = attribute_value;
   }

   eval.checkVariable(var);

   G4int _var  = eval.EvaluateInteger(var);
   G4int _from = eval.EvaluateInteger(from);
   G4int _to   = eval.EvaluateInteger(to);
   G4int _step = eval.EvaluateInteger(step);
   
   if (!from.empty()) _var = _from;

   while (_var <= _to) {
   
      eval.setVariable(var,_var);
      volume_contentRead(element,volumePtr);

      _var += _step;
   }
}

void G4GDMLStructure::structureRead(const xercesc::DOMElement* const element) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="loop") loopRead(child); else
      if (tag=="volume") volumeRead(child); else      
      G4Exception("GDML: Unknown tag in structure: "+tag);
   }
}

G4LogicalVolume* G4GDMLStructure::getVolume(const G4String& ref) const {

   G4LogicalVolume *volumePtr = G4LogicalVolumeStore::GetInstance()->GetVolume(ref,false);

   if (!volumePtr) G4Exception("GDML: Referenced volume '"+ref+"' was not found!");

   return volumePtr;
}
