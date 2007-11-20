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
// $Id: G4GDMLStructure.cc,v 1.11 2007-11-20 13:54:05 ztorzsok Exp $
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

G4bool G4GDMLStructure::directionRead(const xercesc::DOMElement* const element,EAxis& axis) {

   G4String x,y,z;

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

   G4double _x,_y,_z;

   if (!evaluator->Evaluate(_x,x)) return false;
   if (!evaluator->Evaluate(_y,y)) return false;
   if (!evaluator->Evaluate(_z,z)) return false;

   if (_x == 1.0 && _y == 0.0 && _z == 0.0) { axis = kXAxis; return true; } else
   if (_x == 0.0 && _y == 1.0 && _z == 0.0) { axis = kYAxis; return true; } else
   if (_x == 0.0 && _y == 0.0 && _z == 1.0) { axis = kZAxis; return true; }

   G4cout << "GDML ERROR! Only directions along axes are supported!"  << G4endl;

   return false;
}

G4bool G4GDMLStructure::divisionvolRead(const xercesc::DOMElement* const element,G4LogicalVolume* pMother) {

   G4String unit,axis,width,offset,number,volumeref;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="unit"  ) { unit   = attribute_value; } else
      if (attribute_name=="axis"  ) { unit   = attribute_value; } else
      if (attribute_name=="width" ) { width  = attribute_value; } else
      if (attribute_name=="offset") { offset = attribute_value; } else
      if (attribute_name=="number") { number = attribute_value; }
   }

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="volumeref") { if (!refRead(child,volumeref)) return false; } 
   }

   G4LogicalVolume* pLogical = getVolume(module+volumeref);

   if (!pLogical) return false;

   EAxis    _axis   = kZAxis;
   G4double _width  = 0.0;
   G4double _offset = 0.0;
   G4double _number = 0.0;
   
   if (axis=="kXAxis") { _axis = kXAxis; } else
   if (axis=="kYAxis") { _axis = kYAxis; } else
   if (axis=="kZAxis") { _axis = kZAxis; } else
   if (axis=="kRho"  ) { _axis = kRho;   } else
   if (axis=="kPhi"  ) { _axis = kPhi;   }

   if (!evaluator->Evaluate(_width ,width ,unit)) return false;
   if (!evaluator->Evaluate(_offset,offset,unit)) return false;
   if (!evaluator->Evaluate(_number,number))      return false;

   new G4PVDivision("",pLogical,pMother,_axis,(G4int)_number,_width,_offset);

   return true;
}

G4bool G4GDMLStructure::fileRead(const xercesc::DOMElement* const element,G4LogicalVolume** logvol) {

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

      if (attribute_name=="name"   ) { name    = attribute_value; } else
      if (attribute_name=="volname") { volname = attribute_value; }
   }

   G4GDMLStructure structure;
   
   if (!structure.gdmlRead(name,parser)) return false;

   *logvol = structure.getVolume(structure.module + volname);

   return (*logvol != 0);
}

G4bool G4GDMLStructure::loopRead(const xercesc::DOMElement* const element) {

   G4String ffor,to,step;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="for" ) { ffor = attribute_value; } else
      if (attribute_name=="to"  ) { to   = attribute_value; } else
      if (attribute_name=="step") { step = attribute_value; }
   }

   return true;
}

G4bool G4GDMLStructure::paramvolRead(const xercesc::DOMElement* const element,G4LogicalVolume *pMotherLogical) {

   G4String volumeref,parameterised_position_size;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="parameterised_position_size") { ; }  else
      if (tag=="volumeref"                  ) { if (!refRead(child,volumeref)) return false; } 
   }

   return true;
}

G4bool G4GDMLStructure::physvolRead(const xercesc::DOMElement* const element,G4LogicalVolume *mother) {

   G4String volumeref,positionref,rotationref;

   G4LogicalVolume* logvol = 0;
   G4ThreeVector tlate;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="file"       ) { if (!fileRead    (child,&logvol    )) return false; } else
      if (tag=="volumeref"  ) { if (!refRead     (child,volumeref  )) return false; } else
      if (tag=="position"   ) { if (!positionRead(child,tlate      )) return false; } else
      if (tag=="positionref") { if (!refRead     (child,positionref)) return false; } else
      if (tag=="rotationref") { if (!refRead     (child,rotationref)) return false; }
   }

   if (!volumeref.empty()) {

      logvol = getVolume(module+volumeref);

      if (!logvol) return false;
   }

   if (!positionref.empty()) {
   
      G4ThreeVector* posPtr = solids.define.getPosition(module+positionref);

      if (!posPtr) return false;

      tlate = *posPtr;            
   }

   G4RotationMatrix* pRot=0;

   if (!rotationref.empty()) {

      G4ThreeVector* anglePtr = solids.define.getRotation(module+rotationref);

      if (!anglePtr) return false;

      pRot = new G4RotationMatrix();

      pRot->rotateX(anglePtr->x());
      pRot->rotateY(anglePtr->y());
      pRot->rotateZ(anglePtr->z());
   }
 
   new G4PVPlacement(pRot,tlate,logvol,"",mother,false,0);
 
   return true;
}

G4bool G4GDMLStructure::positionRead(const xercesc::DOMElement* const element,G4ThreeVector& vect) {

   G4String unit,x,y,z;

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

   G4double _x,_y,_z;

   if (!evaluator->Evaluate(_x,x,unit)) return false;
   if (!evaluator->Evaluate(_y,y,unit)) return false;
   if (!evaluator->Evaluate(_z,z,unit)) return false;

   vect.set(_x,_y,_z);

   return true;
}

G4bool G4GDMLStructure::quantityRead(const xercesc::DOMElement* const element,G4double& _value) {

   G4String value,unit;

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

   return evaluator->Evaluate(_value,value,unit);
}

G4bool G4GDMLStructure::refRead(const xercesc::DOMElement* const element,G4String& ref) {

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="ref") { ref = attribute_value; }
   }

   return true;
}

G4bool G4GDMLStructure::replicate_along_axisRead(const xercesc::DOMElement* const element,G4double& _width,G4double& _offset,EAxis& _axis) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="width"    ) { if (!quantityRead(child,_width))  return false; } else
      if (tag=="offset"   ) { if (!quantityRead(child,_offset)) return false; } else
      if (tag=="direction") { if (!directionRead(child,_axis))  return false; }
   }

   return true;
}

G4bool G4GDMLStructure::replicavolRead(const xercesc::DOMElement* const element,G4LogicalVolume* pMother) {

   G4String volumeref,numb;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="numb") { numb = attribute_value; }
   }

   G4double _numb;
   G4double _width = 0.0;
   G4double _offset = 0.0;
   EAxis _axis;

   if (!evaluator->Evaluate(_numb,numb)) return false;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="volumeref"           ) { if (!refRead(child,volumeref))                             return false; } else
      if (tag=="replicate_along_axis") { if (!replicate_along_axisRead(child,_width,_offset,_axis)) return false; }
   }

   G4LogicalVolume* pLogical = getVolume(module+volumeref);

   if (!pLogical) return false;

   new G4PVReplica("",pLogical,pMother,_axis,(G4int)_numb,_width,_offset);

   return true;
}

G4bool G4GDMLStructure::volumeRead(const xercesc::DOMElement* const element) {

   G4String name,solidref,materialref;

   XMLCh *name_attr = xercesc::XMLString::transcode("name");
   name = xercesc::XMLString::transcode(element->getAttribute(name_attr));
   xercesc::XMLString::release(&name_attr);

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="materialref") { if (!refRead(child,materialref)) return false; } else
      if (tag=="solidref"   ) { if (!refRead(child,solidref   )) return false; }
   }

   G4Material *materialPtr = materials.getMaterial(module+materialref);
   G4VSolid* solidPtr = solids.getSolid(module+solidref);

   if (!solidPtr || !materialPtr) return false;
 
   G4LogicalVolume *volumePtr = new G4LogicalVolume(solidPtr,materialPtr,module+name,0,0,0);
   
   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);
  
      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="paramvol"   ) { if (!paramvolRead   (child,volumePtr)) return false; } else
      if (tag=="physvol"    ) { if (!physvolRead    (child,volumePtr)) return false; } else
      if (tag=="replicavol" ) { if (!replicavolRead (child,volumePtr)) return false; } else
      if (tag=="divisionvol") { if (!divisionvolRead(child,volumePtr)) return false; }
   }

   return true;
}

G4bool G4GDMLStructure::Read(const xercesc::DOMElement* const element) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="loop"  ) { if (!loopRead  (child)) return false; } else
      if (tag=="volume") { if (!volumeRead(child)) return false; } else      
      {
	 G4cout << "GDML: Error! Unknown tag in structure: " << tag << G4endl;
         return false;
      }
   }
  
   return true;
}

G4bool G4GDMLStructure::gdmlRead(const G4String& fileName,xercesc::XercesDOMParser* newParser) {

   module = fileName + "::";

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
   catch (...) {

      G4cout << "Unexpected exception!" << G4endl;
   }

   xercesc::DOMDocument* doc = parser->getDocument();

   if (!doc) {
   
      G4cout << "Unable to open document '" << fileName << "'!" << G4endl;
      return false;
   }

   xercesc::DOMElement* element = doc->getDocumentElement();

   if (!element) {
   
      G4cout << "Empty document!" << G4endl;
      return false;
   }

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="define"   ) { if (!solids.define.Read(child,evaluator,module)) return false; } else
      if (tag=="materials") { if (!materials.Read    (child,evaluator,module)) return false; } else
      if (tag=="solids"   ) { if (!solids.Read       (child,evaluator,module,"")) return false; } else
      if (tag=="setup"    ) { if (!setup.Read        (child,module)) return false; } else
      if (tag=="structure") { if (!Read              (child)) return false; }
   }

   return true;
}

G4LogicalVolume* G4GDMLStructure::getVolume(const G4String& ref) const {

   G4LogicalVolume *volumePtr = G4LogicalVolumeStore::GetInstance()->GetVolume(ref,false);

   if (!volumePtr) G4cout << "GDML: Error! Referenced volume '" << ref << "' was not found!" << G4endl;   

   return volumePtr;
}
