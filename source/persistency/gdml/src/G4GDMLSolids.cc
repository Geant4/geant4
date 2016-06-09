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
// $Id: G4GDMLSolids.cc,v 1.23.2.1 2008/01/16 09:43:09 ztorzsok Exp $
// GEANT4 tag $ Name:$
//
// class G4GDMLSolids Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLSolids.hh"

void G4GDMLSolids::booleanRead(const xercesc::DOMElement* const element,const BooleanOp op) {

   G4String name;
   G4String first;
   G4String second;

   G4ThreeVector position;
   G4ThreeVector rotation;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value;
   }

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="first") first = refRead(child); else
      if (tag=="second") second = refRead(child); else
      if (tag=="position") position = positionRead(child); else
      if (tag=="rotation") rotation = rotationRead(child); else
      G4Exception("GDML: Unknown tag in boolean solid: "+tag);
   }

   G4RotationMatrix rot;

   rot.rotateX(rotation.x());
   rot.rotateY(rotation.y());
   rot.rotateZ(rotation.z());

   G4Transform3D transform(rot,position);

   G4VSolid* firstSolid = getSolid(GenerateName(first));
   G4VSolid* secondSolid = getSolid(GenerateName(second));

   if (op==UNION) new G4UnionSolid(GenerateName(name),firstSolid,secondSolid,transform); else
   if (op==SUBTRACTION) new G4SubtractionSolid(GenerateName(name),firstSolid,secondSolid,transform); else
   if (op==INTERSECTION) new G4IntersectionSolid(GenerateName(name),firstSolid,secondSolid,transform);
}

void G4GDMLSolids::boxRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit("1");
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

      if (attribute_name=="name" ) name  = attribute_value; else
      if (attribute_name=="lunit") lunit = attribute_value; else
      if (attribute_name=="x") x = attribute_value; else
      if (attribute_name=="y") y = attribute_value; else
      if (attribute_name=="z") z = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);

   G4double _x = eval.Evaluate(x)*_lunit;
   G4double _y = eval.Evaluate(y)*_lunit;
   G4double _z = eval.Evaluate(z)*_lunit;

   _x *= 0.5;
   _y *= 0.5;
   _z *= 0.5;

   new G4Box(GenerateName(name),_x,_y,_z);
}

void G4GDMLSolids::coneRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit("1");
   G4String aunit("1");
   G4String rmin1;
   G4String rmax1;
   G4String rmin2;
   G4String rmax2;
   G4String z;
   G4String startphi;
   G4String deltaphi;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="lunit") lunit = attribute_value; else
      if (attribute_name=="aunit") aunit = attribute_value; else
      if (attribute_name=="rmin1") rmin1 = attribute_value; else
      if (attribute_name=="rmax1") rmax1 = attribute_value; else
      if (attribute_name=="rmin2") rmin2 = attribute_value; else
      if (attribute_name=="rmax2") rmax2 = attribute_value; else
      if (attribute_name=="z") z = attribute_value; else
      if (attribute_name=="startphi") startphi = attribute_value; else
      if (attribute_name=="deltaphi") deltaphi = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);
   G4double _aunit = eval.Evaluate(aunit);

   G4double _rmin1 = eval.Evaluate(rmin1)*_lunit;
   G4double _rmax1 = eval.Evaluate(rmax1)*_lunit;
   G4double _rmin2 = eval.Evaluate(rmin2)*_lunit;
   G4double _rmax2 = eval.Evaluate(rmax2)*_lunit;
   G4double _z = eval.Evaluate(z)*_lunit;
   G4double _startphi = eval.Evaluate(startphi)*_aunit;
   G4double _deltaphi = eval.Evaluate(deltaphi)*_aunit;;

   _z *= 0.5;

   new G4Cons(GenerateName(name),_rmin1,_rmax1,_rmin2,_rmax2,_z,_startphi,_deltaphi);
}

void G4GDMLSolids::ellipsoidRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit("1");
   G4String ax;
   G4String by;
   G4String cz;
   G4String zcut1;
   G4String zcut2; 

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name  = attribute_value; else
      if (attribute_name=="lunit") lunit = attribute_value; else
      if (attribute_name=="ax") ax = attribute_value; else
      if (attribute_name=="by") by = attribute_value; else
      if (attribute_name=="cz") cz = attribute_value; else
      if (attribute_name=="zcut1") zcut1 = attribute_value; else
      if (attribute_name=="zcut2") zcut2 = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);

   G4double _ax = eval.Evaluate(ax)*_lunit;
   G4double _by = eval.Evaluate(by)*_lunit;
   G4double _cz = eval.Evaluate(cz)*_lunit;
   G4double _zcut1 = eval.Evaluate(zcut1)*_lunit;
   G4double _zcut2 = eval.Evaluate(zcut2)*_lunit; 

   new G4Ellipsoid(GenerateName(name),_ax,_by,_cz,_zcut1,_zcut2);
}

void G4GDMLSolids::eltubeRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit("1");
   G4String dx;
   G4String dy;
   G4String dz;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name  = attribute_value; else
      if (attribute_name=="lunit") lunit = attribute_value; else
      if (attribute_name=="dx") dx = attribute_value; else
      if (attribute_name=="dy") dy = attribute_value; else
      if (attribute_name=="dz") dz = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);

   G4double _dx = eval.Evaluate(dx)*_lunit;
   G4double _dy = eval.Evaluate(dy)*_lunit;
   G4double _dz = eval.Evaluate(dz)*_lunit;

   new G4EllipticalTube(GenerateName(name),_dx,_dy,_dz);
}

void G4GDMLSolids::hypeRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit("1");
   G4String aunit("1");
   G4String rmin;
   G4String rmax;
   G4String inst;
   G4String outst;
   G4String z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name  = attribute_value; else
      if (attribute_name=="lunit") lunit = attribute_value; else
      if (attribute_name=="aunit") aunit = attribute_value; else
      if (attribute_name=="rmin") rmin = attribute_value; else
      if (attribute_name=="rmax") rmax = attribute_value; else
      if (attribute_name=="inst") inst = attribute_value; else
      if (attribute_name=="outst") outst = attribute_value; else
      if (attribute_name=="z") z = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);
   G4double _aunit = eval.Evaluate(aunit);

   G4double _rmin = eval.Evaluate(rmin)*_lunit;
   G4double _rmax = eval.Evaluate(rmax)*_lunit;;
   G4double _inst = eval.Evaluate(inst)*_aunit;
   G4double _outst = eval.Evaluate(outst)*_aunit;
   G4double _z = eval.Evaluate(z)*_lunit;

   _z *= 0.5;

   new G4Hype(GenerateName(name),_rmin,_rmax,_inst,_outst,_z);
}

void G4GDMLSolids::loopRead(const xercesc::DOMElement* const element) {

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

      if (attribute_name=="var") var = attribute_value; else
      if (attribute_name=="from") from = attribute_value; else
      if (attribute_name=="to") to = attribute_value; else
      if (attribute_name=="step") step = attribute_value;
   }

   eval.checkVariable(var);

   G4int _var = eval.EvaluateInteger(var);
   G4int _from = eval.EvaluateInteger(from);
   G4int _to = eval.EvaluateInteger(to);
   G4int _step = eval.EvaluateInteger(step);
   
   if (!from.empty()) _var = _from;

   while (_var <= _to) {
   
      eval.setVariable(var,_var);
      solidsRead(element);

      _var += _step;
   }
}

void G4GDMLSolids::orbRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit("1");
   G4String r;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name  = attribute_value; else
      if (attribute_name=="lunit") lunit = attribute_value; else
      if (attribute_name=="r") r = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);
   
   G4double _r = eval.Evaluate(r)*_lunit;

   new G4Orb(GenerateName(name),_r);
}

void G4GDMLSolids::paraRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit("1");
   G4String aunit("1");
   G4String x;
   G4String y;
   G4String z;
   G4String alpha;
   G4String theta;
   G4String phi;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name" ) name  = attribute_value; else
      if (attribute_name=="lunit") lunit = attribute_value; else
      if (attribute_name=="aunit") aunit = attribute_value; else
      if (attribute_name=="x") x = attribute_value; else
      if (attribute_name=="y") y = attribute_value; else
      if (attribute_name=="z") z = attribute_value; else
      if (attribute_name=="alpha") alpha = attribute_value; else
      if (attribute_name=="theta") theta = attribute_value; else
      if (attribute_name=="phi") phi = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);
   G4double _aunit = eval.Evaluate(aunit);

   G4double _x = eval.Evaluate(x)*_lunit;
   G4double _y = eval.Evaluate(y)*_lunit;
   G4double _z = eval.Evaluate(z)*_lunit;
   G4double _alpha = eval.Evaluate(alpha)*_aunit;
   G4double _theta = eval.Evaluate(theta)*_aunit;
   G4double _phi = eval.Evaluate(phi)*_aunit;

   _x *= 0.5;
   _y *= 0.5;
   _z *= 0.5;

   new G4Para(GenerateName(name),_x,_y,_z,_alpha,_theta,_phi);
}

void G4GDMLSolids::polyconeRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit("1");
   G4String aunit("1");
   G4String startphi;
   G4String deltaphi;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="lunit") lunit = attribute_value; else
      if (attribute_name=="aunit") aunit = attribute_value; else
      if (attribute_name=="startphi") startphi = attribute_value; else
      if (attribute_name=="deltaphi") deltaphi = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);
   G4double _aunit = eval.Evaluate(aunit);

   G4double _startphi = eval.Evaluate(startphi)*_aunit;
   G4double _deltaphi = eval.Evaluate(deltaphi)*_aunit;

   std::vector<zplaneType> zplaneList;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="zplane") zplaneList.push_back(zplaneRead(child,_lunit));
   }

   G4int numZPlanes = zplaneList.size();

   G4double *rmin_array = new G4double[numZPlanes];
   G4double *rmax_array = new G4double[numZPlanes];
   G4double* z_array    = new G4double[numZPlanes];

   for (G4int i=0;i<numZPlanes;i++) {
   
      rmin_array[i] = zplaneList[i].rmin;
      rmax_array[i] = zplaneList[i].rmax;
      z_array[i]    = zplaneList[i].z;
   }

   new G4Polycone(GenerateName(name),_startphi,_deltaphi,numZPlanes,z_array,rmin_array,rmax_array);
}

void G4GDMLSolids::polyhedraRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit("1");
   G4String aunit("1");
   G4String startphi;
   G4String deltaphi;
   G4String numsides;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name     = attribute_value; else
      if (attribute_name=="lunit") lunit    = attribute_value; else
      if (attribute_name=="aunit" ) aunit    = attribute_value; else
      if (attribute_name=="startphi") startphi = attribute_value; else
      if (attribute_name=="deltaphi") deltaphi = attribute_value; else
      if (attribute_name=="numsides") numsides = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);
   G4double _aunit = eval.Evaluate(aunit);

   G4double _startphi = eval.Evaluate(startphi)*_aunit;
   G4double _deltaphi = eval.Evaluate(deltaphi)*_aunit;

   G4int _numsides = eval.EvaluateInteger(numsides);

   std::vector<zplaneType> zplaneList;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="zplane") zplaneList.push_back(zplaneRead(child,_lunit));
   }

   G4int numZPlanes = zplaneList.size();

   G4double *rmin_array = new G4double[numZPlanes];
   G4double *rmax_array = new G4double[numZPlanes];
   G4double* z_array    = new G4double[numZPlanes];

   for (G4int i=0;i<numZPlanes;i++) {
   
      rmin_array[i] = zplaneList[i].rmin;
      rmax_array[i] = zplaneList[i].rmax;
      z_array[i]    = zplaneList[i].z;
   }

   new G4Polyhedra(GenerateName(name),_startphi,_deltaphi,_numsides,numZPlanes,z_array,rmin_array,rmax_array);
}

G4QuadrangularFacet* G4GDMLSolids::quadrangularRead(const xercesc::DOMElement* const element) {

   G4String v1;
   G4String v2;
   G4String v3;
   G4String v4;
   G4String type;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="vertex1") v1 = attribute_value; else
      if (attribute_name=="vertex2") v2 = attribute_value; else
      if (attribute_name=="vertex3") v3 = attribute_value; else
      if (attribute_name=="vertex4") v4 = attribute_value; else
      if (attribute_name=="type") type = attribute_value;
   }

   const G4ThreeVector* ptr1 = getPosition(GenerateName(v1));
   const G4ThreeVector* ptr2 = getPosition(GenerateName(v2));
   const G4ThreeVector* ptr3 = getPosition(GenerateName(v3));
   const G4ThreeVector* ptr4 = getPosition(GenerateName(v4));

   return new G4QuadrangularFacet(*ptr1,*ptr2,*ptr3,*ptr4,(type=="RELATIVE")?(RELATIVE):(ABSOLUTE));
}

G4String G4GDMLSolids::refRead(const xercesc::DOMElement* const element) {

   G4String ref;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="ref") ref = attribute_value;
   }

   return ref;
}

void G4GDMLSolids::reflectedSolidRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit("1");
   G4String aunit("1");
   G4String solid;
   G4String sx;
   G4String sy;
   G4String sz;
   G4String rx;
   G4String ry;
   G4String rz;
   G4String dx;
   G4String dy;
   G4String dz;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name" ) name  = attribute_value; else
      if (attribute_name=="lunit") lunit = attribute_value; else
      if (attribute_name=="aunit") aunit = attribute_value; else
      if (attribute_name=="solid") solid = attribute_value; else
      if (attribute_name=="sx") sx = attribute_value; else
      if (attribute_name=="sy") sy = attribute_value; else
      if (attribute_name=="sz") sz = attribute_value; else
      if (attribute_name=="rx") rx = attribute_value; else
      if (attribute_name=="ry") ry = attribute_value; else
      if (attribute_name=="rz") rz = attribute_value; else
      if (attribute_name=="dx") dx = attribute_value; else
      if (attribute_name=="dy") dy = attribute_value; else
      if (attribute_name=="dz") dz = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);
   G4double _aunit = eval.Evaluate(aunit);

   G4double _sx = eval.Evaluate(sx);
   G4double _sy = eval.Evaluate(sy);
   G4double _sz = eval.Evaluate(sz);
   G4double _rx = eval.Evaluate(rx)*_aunit;
   G4double _ry = eval.Evaluate(ry)*_aunit;
   G4double _rz = eval.Evaluate(rz)*_aunit;
   G4double _dx = eval.Evaluate(dx)*_lunit;
   G4double _dy = eval.Evaluate(dy)*_lunit;
   G4double _dz = eval.Evaluate(dz)*_lunit;

   G4VSolid* solidPtr = getSolid(GenerateName(solid));

   G4RotationMatrix rot;
   
   rot.rotateX(_rx);
   rot.rotateY(_ry);
   rot.rotateZ(_rz);

   G4ThreeVector trans(_dx,_dy,_dz);

   G4Scale3D scale(_sx,_sy,_sz);

   G4Transform3D transform(rot,trans);
   transform = transform*scale;
          
   new G4ReflectedSolid(GenerateName(name),solidPtr,transform);
}

G4ExtrudedSolid::ZSection G4GDMLSolids::sectionRead(const xercesc::DOMElement* const element,G4double _lunit) {

   G4String zPosition;
   G4String xOffset;
   G4String yOffset;
   G4String scalingFactor;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="zPosition") zPosition = attribute_value; else
      if (attribute_name=="xOffset") xOffset = attribute_value; else
      if (attribute_name=="yOffset") yOffset = attribute_value; else
      if (attribute_name=="scalingFactor") scalingFactor = attribute_value;
   }

   G4double _zPosition = eval.Evaluate(zPosition)*_lunit;
   G4double _xOffset = eval.Evaluate(xOffset)*_lunit;
   G4double _yOffset = eval.Evaluate(yOffset)*_lunit;
   G4double _scalingFactor = eval.Evaluate(scalingFactor);

   return G4ExtrudedSolid::ZSection(_zPosition,G4TwoVector(_xOffset,_yOffset),_scalingFactor);
}

void G4GDMLSolids::sphereRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit("1");
   G4String aunit("1");
   G4String rmin;
   G4String rmax;
   G4String startphi;
   G4String deltaphi;
   G4String starttheta;
   G4String deltatheta;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="lunit") lunit = attribute_value; else
      if (attribute_name=="aunit") aunit = attribute_value; else
      if (attribute_name=="rmin") rmin = attribute_value; else
      if (attribute_name=="rmax") rmax = attribute_value; else
      if (attribute_name=="startphi") startphi = attribute_value; else
      if (attribute_name=="deltaphi") deltaphi = attribute_value; else
      if (attribute_name=="starttheta") starttheta = attribute_value; else
      if (attribute_name=="deltatheta") deltatheta = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);
   G4double _aunit = eval.Evaluate(aunit);

   G4double _rmin = eval.Evaluate(rmin)*_lunit;
   G4double _rmax = eval.Evaluate(rmax)*_lunit;
   G4double _startphi = eval.Evaluate(startphi)*_aunit;
   G4double _deltaphi = eval.Evaluate(deltaphi)*_aunit;
   G4double _starttheta = eval.Evaluate(starttheta)*_aunit;
   G4double _deltatheta = eval.Evaluate(deltatheta)*_aunit;

   new G4Sphere(GenerateName(name),_rmin,_rmax,_startphi,_deltaphi,_starttheta,_deltatheta);
}

void G4GDMLSolids::tessellatedRead(const xercesc::DOMElement* const element) {

   G4String name;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value;
   }
   
   G4TessellatedSolid *tessellated = new G4TessellatedSolid(GenerateName(name));

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="triangular") tessellated->AddFacet(triangularRead(child)); else
      if (tag=="quadrangular") tessellated->AddFacet(quadrangularRead(child));
   }

   tessellated->SetSolidClosed(true);
}

void G4GDMLSolids::tetRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String vertex1;
   G4String vertex2;
   G4String vertex3;
   G4String vertex4;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="vertex1") vertex1 = attribute_value; else
      if (attribute_name=="vertex2") vertex2 = attribute_value; else
      if (attribute_name=="vertex3") vertex3 = attribute_value; else
      if (attribute_name=="vertex4") vertex4 = attribute_value;
   }
   
   const G4ThreeVector* ptr1 = getPosition(GenerateName(vertex1));
   const G4ThreeVector* ptr2 = getPosition(GenerateName(vertex2));
   const G4ThreeVector* ptr3 = getPosition(GenerateName(vertex3));
   const G4ThreeVector* ptr4 = getPosition(GenerateName(vertex4));

   new G4Tet(GenerateName(name),*ptr1,*ptr2,*ptr3,*ptr4);
}

void G4GDMLSolids::torusRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit("1");
   G4String aunit("1");
   G4String rmin;
   G4String rmax;
   G4String rtor;
   G4String startphi;
   G4String deltaphi;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="lunit") lunit = attribute_value; else
      if (attribute_name=="aunit") aunit = attribute_value; else
      if (attribute_name=="rmin") rmin = attribute_value; else
      if (attribute_name=="rmax") rmax = attribute_value; else
      if (attribute_name=="rtor") rtor = attribute_value; else
      if (attribute_name=="startphi") startphi = attribute_value; else
      if (attribute_name=="deltaphi") deltaphi = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);
   G4double _aunit = eval.Evaluate(aunit);

   G4double _rmin = eval.Evaluate(rmin)*_lunit;
   G4double _rmax = eval.Evaluate(rmax)*_lunit;
   G4double _rtor = eval.Evaluate(rtor)*_lunit;
   G4double _startphi = eval.Evaluate(startphi)*_aunit;
   G4double _deltaphi = eval.Evaluate(deltaphi)*_aunit;

   new G4Torus(GenerateName(name),_rmin,_rmax,_rtor,_startphi,_deltaphi);
}

void G4GDMLSolids::trapRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit("1");
   G4String aunit("1");
   G4String z;
   G4String theta;
   G4String phi;
   G4String y1;
   G4String x1;
   G4String x2;
   G4String alpha1;
   G4String y2;
   G4String x3;
   G4String x4;
   G4String alpha2;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="lunit") lunit = attribute_value; else
      if (attribute_name=="aunit") aunit = attribute_value; else
      if (attribute_name=="z") z = attribute_value; else
      if (attribute_name=="theta") theta  = attribute_value; else
      if (attribute_name=="phi") phi = attribute_value; else
      if (attribute_name=="y1") y1 = attribute_value; else
      if (attribute_name=="x1") x1 = attribute_value; else
      if (attribute_name=="x2") x2 = attribute_value; else
      if (attribute_name=="alpha1") alpha1 = attribute_value; else
      if (attribute_name=="y2") y2 = attribute_value; else
      if (attribute_name=="x3") x3 = attribute_value; else
      if (attribute_name=="x4") x4 = attribute_value; else
      if (attribute_name=="alpha2") alpha2 = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);
   G4double _aunit = eval.Evaluate(aunit);

   G4double _z = eval.Evaluate(z)*_lunit;
   G4double _theta  = eval.Evaluate(theta)*_aunit;
   G4double _phi = eval.Evaluate(phi)*_aunit;
   G4double _y1 = eval.Evaluate(y1)*_lunit;
   G4double _x1 = eval.Evaluate(x1)*_lunit;
   G4double _x2 = eval.Evaluate(x2)*_lunit;
   G4double _alpha1 = eval.Evaluate(alpha1)*_aunit;
   G4double _y2 = eval.Evaluate(y2)*_lunit;
   G4double _x3 = eval.Evaluate(x3)*_lunit;
   G4double _x4 = eval.Evaluate(x4)*_lunit;
   G4double _alpha2 = eval.Evaluate(alpha2)*_aunit;

   _z *= 0.5;
   _y1 *= 0.5;
   _x1 *= 0.5;
   _x2 *= 0.5;
   _y2 *= 0.5;
   _x3 *= 0.5;
   _x4 *= 0.5;

   new G4Trap(GenerateName(name),_z,_theta,_phi,_y1,_x1,_x2,_alpha1,_y2,_x3,_x4,_alpha2);
}

void G4GDMLSolids::trdRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit("1");
   G4String x1;
   G4String x2;
   G4String y1;
   G4String y2;
   G4String z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="lunit") lunit = attribute_value; else
      if (attribute_name=="x1") x1 = attribute_value; else
      if (attribute_name=="x2") x2 = attribute_value; else
      if (attribute_name=="y1") y1 = attribute_value; else
      if (attribute_name=="y2") y2 = attribute_value; else
      if (attribute_name=="z") z = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);

   G4double _x1 = eval.Evaluate(x1)*_lunit;
   G4double _x2 = eval.Evaluate(x2)*_lunit;
   G4double _y1 = eval.Evaluate(y1)*_lunit;
   G4double _y2 = eval.Evaluate(y2)*_lunit;
   G4double _z  = eval.Evaluate(z)*_lunit;

   _x1 *= 0.5;
   _x2 *= 0.5;
   _y1 *= 0.5;
   _y2 *= 0.5;
   _z *= 0.5;

   new G4Trd(GenerateName(name),_x1,_x2,_y1,_y2,_z);
}

G4TriangularFacet* G4GDMLSolids::triangularRead(const xercesc::DOMElement* const element) {

   G4String v1;
   G4String v2;
   G4String v3;
   G4String type;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="vertex1") v1 = attribute_value; else
      if (attribute_name=="vertex2") v2 = attribute_value; else
      if (attribute_name=="vertex3") v3 = attribute_value; else
      if (attribute_name=="type") type = attribute_value;
   }

   const G4ThreeVector* ptr1 = getPosition(GenerateName(v1));
   const G4ThreeVector* ptr2 = getPosition(GenerateName(v2));
   const G4ThreeVector* ptr3 = getPosition(GenerateName(v3));

   return new G4TriangularFacet(*ptr1,*ptr2,*ptr3,(type=="RELATIVE")?(RELATIVE):(ABSOLUTE));
}

void G4GDMLSolids::tubeRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit("1");
   G4String aunit("1");
   G4String rmin;
   G4String rmax;
   G4String z;
   G4String startphi;
   G4String deltaphi;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="lunit") lunit = attribute_value; else
      if (attribute_name=="aunit") aunit = attribute_value; else
      if (attribute_name=="rmin") rmin = attribute_value; else
      if (attribute_name=="rmax") rmax = attribute_value; else
      if (attribute_name=="z") z = attribute_value; else
      if (attribute_name=="startphi") startphi = attribute_value; else
      if (attribute_name=="deltaphi") deltaphi = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);
   G4double _aunit = eval.Evaluate(aunit);

   G4double _rmin = eval.Evaluate(rmin)*_lunit;
   G4double _rmax = eval.Evaluate(rmax)*_lunit;
   G4double _z = eval.Evaluate(z)*_lunit;
   G4double _startphi = eval.Evaluate(startphi)*_aunit;
   G4double _deltaphi = eval.Evaluate(deltaphi)*_aunit;

   _z *= 0.5;

   new G4Tubs(GenerateName(name),_rmin,_rmax,_z,_startphi,_deltaphi);
}

G4TwoVector G4GDMLSolids::twoDimVertexRead(const xercesc::DOMElement* const element,G4double _lunit) {

   G4String x;
   G4String y;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="x") x = attribute_value; else
      if (attribute_name=="y") y = attribute_value;
   }

   G4double _x = eval.Evaluate(x)*_lunit;
   G4double _y = eval.Evaluate(y)*_lunit;

   return G4TwoVector(_x,_y);
}

void G4GDMLSolids::xtruRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit("1");

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") name = attribute_value; else
      if (attribute_name=="lunit") lunit = attribute_value;
   }

   G4double _lunit = eval.Evaluate(lunit);

   std::vector<G4TwoVector> twoDimVertexList;
   std::vector<G4ExtrudedSolid::ZSection> sectionList;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="twoDimVertex") twoDimVertexList.push_back(twoDimVertexRead(child,_lunit)); else
      if (tag=="section") sectionList.push_back(sectionRead(child,_lunit));      
   }

   new G4ExtrudedSolid(GenerateName(name),twoDimVertexList,sectionList);
}

G4GDMLSolids::zplaneType G4GDMLSolids::zplaneRead(const xercesc::DOMElement* const element,G4double _lunit) {

   G4String rmin;
   G4String rmax;
   G4String z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="rmin") rmin = attribute_value; else
      if (attribute_name=="rmax") rmax = attribute_value; else
      if (attribute_name=="z") z = attribute_value;
   }

   zplaneType zplane;

   zplane.rmin = eval.Evaluate(rmin)*_lunit;
   zplane.rmax = eval.Evaluate(rmax)*_lunit;
   zplane.z = eval.Evaluate(z)*_lunit;

   return zplane;
}

void G4GDMLSolids::solidsRead(const xercesc::DOMElement* const element) {

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="box") boxRead(child); else
      if (tag=="cone") coneRead(child); else
      if (tag=="ellipsoid") ellipsoidRead(child); else
      if (tag=="eltube") eltubeRead(child); else
      if (tag=="hype") hypeRead(child); else
      if (tag=="loop") loopRead(child); else
      if (tag=="orb") orbRead(child); else
      if (tag=="para") paraRead(child); else
      if (tag=="polycone") polyconeRead(child); else
      if (tag=="polyhedra") polyhedraRead(child); else
      if (tag=="sphere") sphereRead(child); else
      if (tag=="reflectedSolid") reflectedSolidRead(child); else
      if (tag=="tessellated") tessellatedRead(child); else
      if (tag=="tet") tetRead(child); else
      if (tag=="torus") torusRead(child); else
      if (tag=="trap") trapRead(child); else
      if (tag=="trd") trdRead(child); else
      if (tag=="tube") tubeRead(child); else
      if (tag=="xtru") xtruRead(child); else
      if (tag=="intersection") booleanRead(child,INTERSECTION); else
      if (tag=="subtraction") booleanRead(child,SUBTRACTION); else
      if (tag=="union") booleanRead(child,UNION); else
      G4Exception("GDML: Unknown tag in solids: "+tag);
   }
}

G4ThreeVector G4GDMLSolids::positionRead(const xercesc::DOMElement* const element) {

   G4String unit("1");
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

      if (attribute_name=="unit") unit = attribute_value; else
      if (attribute_name=="x") x  = attribute_value; else
      if (attribute_name=="y") y = attribute_value; else
      if (attribute_name=="z") z = attribute_value;
   }

   G4double _unit = eval.Evaluate(unit);

   G4double _x = eval.Evaluate(x)*_unit;
   G4double _y = eval.Evaluate(y)*_unit;
   G4double _z = eval.Evaluate(z)*_unit;
   
   return G4ThreeVector(_x,_y,_z);
}

G4ThreeVector G4GDMLSolids::rotationRead(const xercesc::DOMElement* const element) {

   G4String unit("1");
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

      if (attribute_name=="unit") unit = attribute_value; else
      if (attribute_name=="x") x = attribute_value; else
      if (attribute_name=="y") y = attribute_value; else
      if (attribute_name=="z") z = attribute_value;
   }

   G4double _unit = eval.Evaluate(unit);

   G4double _x = eval.Evaluate(x)*_unit;
   G4double _y = eval.Evaluate(y)*_unit;
   G4double _z = eval.Evaluate(z)*_unit;
   
   return G4ThreeVector(_x,_y,_z);
}

G4ThreeVector G4GDMLSolids::scaleRead(const xercesc::DOMElement* const element) {

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
   
   return G4ThreeVector(_x,_y,_z);
}

G4VSolid* G4GDMLSolids::getSolid(const G4String& ref) const {

   G4VSolid* solidPtr = G4SolidStore::GetInstance()->GetSolid(ref,false);

   if (!solidPtr) G4Exception("GDML: Referenced solid '"+ref+"' was not found!");

   return solidPtr;
}
