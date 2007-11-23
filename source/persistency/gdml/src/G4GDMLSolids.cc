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
// $Id: G4GDMLSolids.cc,v 1.16 2007-11-23 10:31:55 ztorzsok Exp $
// GEANT4 tag $ Name:$
//
// class G4GDMLSolids Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLSolids.hh"
#include <sstream>

std::string G4GDMLSolids::nameProcess(const std::string& in) {

   std::string out(file);
   
   std::string::size_type open = in.find("]",0);

   if (open != std::string::npos) G4Exception("Bracket mismatch in loop!");

   open = in.find("[",0);

   out.append(in,0,open);
   
   while (open != std::string::npos) {
   
      std::string::size_type close = in.find("]",open);

      if (close == std::string::npos) G4Exception("Bracket mismatch in loop!");
   
      std::string expr = in.substr(open+1,close-open-1);

      double _expr;
   
      evaluator->Evaluate(_expr,expr);
   
      std::stringstream ss;
   
      ss << "[" << _expr << "]";
   
      out.append(ss.str());

      open = in.find("[",close);
   }

   return out;
}

G4bool G4GDMLSolids::booleanRead(const xercesc::DOMElement* const element,const BooleanOp op) {

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

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") { name = attribute_value; }
   }

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="first"   ) { if (!refRead     (child,first   )) return false; } else
      if (tag=="second"  ) { if (!refRead     (child,second  )) return false; } else
      if (tag=="position") { if (!positionRead(child,position)) return false; } else
      if (tag=="rotation") { if (!rotationRead(child,rotation)) return false; } else
      G4Exception("GDML: Unknown tag in boolean solid: "+tag);
   }

   G4RotationMatrix rot;

   rot.rotateX(rotation.x());
   rot.rotateY(rotation.y());
   rot.rotateZ(rotation.z());

   G4Transform3D transform(rot,position);

   G4VSolid* firstSolid = getSolid(nameProcess(first));
   G4VSolid* secondSolid = getSolid(nameProcess(second));

   if (op==UNION       ) { new G4UnionSolid       (nameProcess(name),firstSolid,secondSolid,transform); } else
   if (op==SUBTRACTION ) { new G4SubtractionSolid (nameProcess(name),firstSolid,secondSolid,transform); } else
   if (op==INTERSECTION) { new G4IntersectionSolid(nameProcess(name),firstSolid,secondSolid,transform); }

   return true;
}

G4bool G4GDMLSolids::boxRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit;
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

      if (attribute_name=="name" ) { name  = attribute_value; } else
      if (attribute_name=="lunit") { lunit = attribute_value; } else
      if (attribute_name=="x"    ) { x     = attribute_value; } else
      if (attribute_name=="y"    ) { y     = attribute_value; } else
      if (attribute_name=="z"    ) { z     = attribute_value; }
   }

   G4double _x;
   G4double _y;
   G4double _z;

   if (!evaluator->Evaluate(_x,x,lunit)) return false;
   if (!evaluator->Evaluate(_y,y,lunit)) return false;
   if (!evaluator->Evaluate(_z,z,lunit)) return false;

   _x *= 0.5;
   _y *= 0.5;
   _z *= 0.5;

   new G4Box(nameProcess(name),_x,_y,_z);

   return true;
}

G4bool G4GDMLSolids::coneRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit;
   G4String aunit;
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

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name"    ) { name     = attribute_value; } else
      if (attribute_name=="lunit"   ) { lunit    = attribute_value; } else
      if (attribute_name=="aunit"   ) { aunit    = attribute_value; } else
      if (attribute_name=="rmin1"   ) { rmin1    = attribute_value; } else
      if (attribute_name=="rmax1"   ) { rmax1    = attribute_value; } else
      if (attribute_name=="rmin2"   ) { rmin2    = attribute_value; } else
      if (attribute_name=="rmax2"   ) { rmax2    = attribute_value; } else
      if (attribute_name=="z"       ) { z        = attribute_value; } else
      if (attribute_name=="startphi") { startphi = attribute_value; } else
      if (attribute_name=="deltaphi") { deltaphi = attribute_value; }
   }

   G4double _rmin1;
   G4double _rmax1;
   G4double _rmin2;
   G4double _rmax2;
   G4double _z;
   G4double _startphi;
   G4double _deltaphi;

   if (!evaluator->Evaluate(_rmin1   ,rmin1   ,lunit)) return false;
   if (!evaluator->Evaluate(_rmax1   ,rmax1   ,lunit)) return false;
   if (!evaluator->Evaluate(_rmin2   ,rmin2   ,lunit)) return false;
   if (!evaluator->Evaluate(_rmax2   ,rmax2   ,lunit)) return false;
   if (!evaluator->Evaluate(_z       ,z       ,lunit)) return false;
   if (!evaluator->Evaluate(_startphi,startphi,aunit)) return false;
   if (!evaluator->Evaluate(_deltaphi,deltaphi,aunit)) return false;

   _z *= 0.5;

   new G4Cons(nameProcess(name),_rmin1,_rmax1,_rmin2,_rmax2,_z,_startphi,_deltaphi);

   return true;
}

G4bool G4GDMLSolids::ellipsoidRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit;
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

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name" ) { name  = attribute_value; } else
      if (attribute_name=="lunit") { lunit = attribute_value; } else
      if (attribute_name=="ax"   ) { ax    = attribute_value; } else
      if (attribute_name=="by"   ) { by    = attribute_value; } else
      if (attribute_name=="cz"   ) { cz    = attribute_value; } else
      if (attribute_name=="zcut1") { zcut1 = attribute_value; } else
      if (attribute_name=="zcut2") { zcut2 = attribute_value; }
   }

   G4double _ax;
   G4double _by;
   G4double _cz;
   G4double _zcut1;
   G4double _zcut2; 

   if (!evaluator->Evaluate(_ax   ,ax   ,lunit)) return false;
   if (!evaluator->Evaluate(_by   ,by   ,lunit)) return false;
   if (!evaluator->Evaluate(_cz   ,cz   ,lunit)) return false;
   if (!evaluator->Evaluate(_zcut1,zcut1,lunit)) return false;
   if (!evaluator->Evaluate(_zcut2,zcut2,lunit)) return false;

   new G4Ellipsoid(nameProcess(name),_ax,_by,_cz,_zcut1,_zcut2);

   return true;
}

G4bool G4GDMLSolids::eltubeRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit;
   G4String dx;
   G4String dy;
   G4String dz;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name" ) { name  = attribute_value; } else
      if (attribute_name=="lunit") { lunit = attribute_value; } else
      if (attribute_name=="dx"   ) { dx    = attribute_value; } else
      if (attribute_name=="dy"   ) { dy    = attribute_value; } else
      if (attribute_name=="dz"   ) { dz    = attribute_value; }
   }

   G4double _dx;
   G4double _dy;
   G4double _dz;

   if (!evaluator->Evaluate(_dx,dx,lunit)) return false;
   if (!evaluator->Evaluate(_dy,dy,lunit)) return false;
   if (!evaluator->Evaluate(_dz,dz,lunit)) return false;

   new G4EllipticalTube(nameProcess(name),_dx,_dy,_dz);

   return true;
}

G4bool G4GDMLSolids::hypeRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit;
   G4String aunit;
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

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name" ) { name  = attribute_value; } else
      if (attribute_name=="lunit") { lunit = attribute_value; } else
      if (attribute_name=="aunit") { aunit = attribute_value; } else
      if (attribute_name=="rmin" ) { rmin  = attribute_value; } else
      if (attribute_name=="rmax" ) { rmax  = attribute_value; } else
      if (attribute_name=="inst" ) { inst  = attribute_value; } else
      if (attribute_name=="outst") { outst = attribute_value; } else
      if (attribute_name=="z"    ) { z     = attribute_value; }
   }

   G4double _rmin;
   G4double _rmax;
   G4double _inst;
   G4double _outst;
   G4double _z;

   if (!evaluator->Evaluate(_rmin ,rmin ,lunit)) return false;
   if (!evaluator->Evaluate(_rmax ,rmax ,lunit)) return false;
   if (!evaluator->Evaluate(_inst ,inst ,aunit)) return false;
   if (!evaluator->Evaluate(_outst,outst,aunit)) return false;
   if (!evaluator->Evaluate(_z    ,z    ,lunit)) return false;

   _z *= 0.5;

   new G4Hype(nameProcess(name),_rmin,_rmax,_inst,_outst,_z);

   return true;
}

G4bool G4GDMLSolids::loopRead(const xercesc::DOMElement* const element) {

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

      if (attribute_name=="var" ) { var  = attribute_value; } else
      if (attribute_name=="from") { from = attribute_value; } else
      if (attribute_name=="to"  ) { to   = attribute_value; } else
      if (attribute_name=="step") { step = attribute_value; }
   }

   G4double _var;
   G4double _from;
   G4double _to;
   G4double _step;
   
   if (!evaluator->Evaluate(_var ,var )) return false;
   if (!evaluator->Evaluate(_from,from)) return false;
   if (!evaluator->Evaluate(_to  ,to  )) return false;
   if (!evaluator->Evaluate(_step,step)) return false;

   if (!from.empty()) _var = _from;

   while (_var <= _to) {
   
      evaluator->setVariable(var,_var);
      Read(element,evaluator,file);

      _var += _step;
   }

   return true;
}

G4bool G4GDMLSolids::orbRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit;
   G4String r;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name" ) { name  = attribute_value; } else
      if (attribute_name=="lunit") { lunit = attribute_value; } else
      if (attribute_name=="r"    ) { r     = attribute_value; }
   }

   G4double _r;

   if (!evaluator->Evaluate(_r,r,lunit)) return false;

   new G4Orb(nameProcess(name),_r);

   return true;
}

G4bool G4GDMLSolids::paraRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit;
   G4String aunit;
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

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name" ) { name  = attribute_value; } else
      if (attribute_name=="lunit") { lunit = attribute_value; } else
      if (attribute_name=="aunit") { aunit = attribute_value; } else
      if (attribute_name=="x"    ) { x     = attribute_value; } else
      if (attribute_name=="y"    ) { y     = attribute_value; } else
      if (attribute_name=="z"    ) { z     = attribute_value; } else
      if (attribute_name=="alpha") { alpha = attribute_value; } else
      if (attribute_name=="theta") { theta = attribute_value; } else
      if (attribute_name=="phi"  ) { phi   = attribute_value; }
   }

   G4double _x;
   G4double _y;
   G4double _z;
   G4double _alpha;
   G4double _theta;
   G4double _phi;

   if (!evaluator->Evaluate(_x    ,x    ,lunit)) return false;
   if (!evaluator->Evaluate(_y    ,y    ,lunit)) return false;
   if (!evaluator->Evaluate(_z    ,z    ,lunit)) return false;
   if (!evaluator->Evaluate(_alpha,alpha,aunit)) return false;
   if (!evaluator->Evaluate(_theta,theta,aunit)) return false;
   if (!evaluator->Evaluate(_phi  ,phi  ,aunit)) return false;

   _x *= 0.5;
   _y *= 0.5;
   _z *= 0.5;

   new G4Para(nameProcess(name),_x,_y,_z,_alpha,_theta,_phi);

   return true;
}

G4bool G4GDMLSolids::polyconeRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit;
   G4String aunit;
   G4String startphi;
   G4String deltaphi;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name"    ) { name     = attribute_value; } else
      if (attribute_name=="lunit"   ) { lunit    = attribute_value; } else
      if (attribute_name=="aunit"   ) { aunit    = attribute_value; } else
      if (attribute_name=="startphi") { startphi = attribute_value; } else
      if (attribute_name=="deltaphi") { deltaphi = attribute_value; }
   }

   G4double _startphi;
   G4double _deltaphi;

   if (!evaluator->Evaluate(_startphi,startphi,aunit)) return false;
   if (!evaluator->Evaluate(_deltaphi,deltaphi,aunit)) return false;

   std::vector<zplaneType> zplaneList;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag!="zplane") continue; 

      zplaneType zplane;      
      
      if (!zplaneRead(child,zplane,lunit)) return false;

      zplaneList.push_back(zplane);
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

   new G4Polycone(nameProcess(name),_startphi,_deltaphi,numZPlanes,z_array,rmin_array,rmax_array);

   return true;
}

G4bool G4GDMLSolids::polyhedraRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit;
   G4String aunit;
   G4String startphi;
   G4String deltaphi;
   G4String numsides;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name"    ) { name     = attribute_value; } else
      if (attribute_name=="lunit"   ) { lunit    = attribute_value; } else
      if (attribute_name=="aunit"   ) { aunit    = attribute_value; } else
      if (attribute_name=="startphi") { startphi = attribute_value; } else
      if (attribute_name=="deltaphi") { deltaphi = attribute_value; } else
      if (attribute_name=="numsides") { numsides = attribute_value; }
   }

   G4double _startphi;
   G4double _deltaphi;
   G4double _numsides;

   if (!evaluator->Evaluate(_startphi,startphi,aunit)) return false;
   if (!evaluator->Evaluate(_deltaphi,deltaphi,aunit)) return false;

   if (!evaluator->Evaluate(_numsides,numsides)) return false;

   std::vector<zplaneType> zplaneList;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="zplane") {

         zplaneType zplane;      
      
         if (!zplaneRead(child,zplane,lunit)) return false;

         zplaneList.push_back(zplane);
      }
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

   new G4Polyhedra(nameProcess(name),_startphi,_deltaphi,(G4int)_numsides,numZPlanes,z_array,rmin_array,rmax_array);
   
   return true;
}

G4bool G4GDMLSolids::positionRead(const xercesc::DOMElement* const element,G4ThreeVector& vect) {

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

   G4double _x;
   G4double _y;
   G4double _z;

   if (!evaluator->Evaluate(_x,x,unit)) return false;
   if (!evaluator->Evaluate(_y,y,unit)) return false;
   if (!evaluator->Evaluate(_z,z,unit)) return false;

   vect.set(_x,_y,_z);

   return true;
}

G4bool G4GDMLSolids::quadrangularRead(const xercesc::DOMElement* const element,G4TessellatedSolid* tessellated) {

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

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="v1"  ) { v1   = attribute_value; } else
      if (attribute_name=="v2"  ) { v2   = attribute_value; } else
      if (attribute_name=="v3"  ) { v3   = attribute_value; } else
      if (attribute_name=="v4"  ) { v4   = attribute_value; } else
      if (attribute_name=="type") { type = attribute_value; }
   }

   const G4ThreeVector* ptr1 = define.getPosition(nameProcess(v1));
   const G4ThreeVector* ptr2 = define.getPosition(nameProcess(v2));
   const G4ThreeVector* ptr3 = define.getPosition(nameProcess(v3));
   const G4ThreeVector* ptr4 = define.getPosition(nameProcess(v4));

   if (!ptr1 || !ptr2 || !ptr3 || !ptr4) return false;

   tessellated->AddFacet(new G4QuadrangularFacet(*ptr1,*ptr2,*ptr3,*ptr4,(type=="RELATIVE")?(RELATIVE):(ABSOLUTE)));

   return true;
}

G4bool G4GDMLSolids::refRead(const xercesc::DOMElement* const element,G4String& ref) {

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

G4bool G4GDMLSolids::reflectedSolidRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit;
   G4String aunit;
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

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name" ) { name  = attribute_value; } else
      if (attribute_name=="lunit") { lunit = attribute_value; } else
      if (attribute_name=="aunit") { aunit = attribute_value; } else
      if (attribute_name=="solid") { solid = attribute_value; } else
      if (attribute_name=="sx"   ) { sx    = attribute_value; } else
      if (attribute_name=="sy"   ) { sy    = attribute_value; } else
      if (attribute_name=="sz"   ) { sz    = attribute_value; } else
      if (attribute_name=="rx"   ) { rx    = attribute_value; } else
      if (attribute_name=="ry"   ) { ry    = attribute_value; } else
      if (attribute_name=="rz"   ) { rz    = attribute_value; } else
      if (attribute_name=="dx"   ) { dx    = attribute_value; } else
      if (attribute_name=="dy"   ) { dy    = attribute_value; } else
      if (attribute_name=="dz"   ) { dz    = attribute_value; }
   }

   G4double _sx;
   G4double _sy;
   G4double _sz;
   G4double _rx;
   G4double _ry;
   G4double _rz;
   G4double _dx;
   G4double _dy;
   G4double _dz;

   if (!evaluator->Evaluate(_sx,sx)) return false; // Scaling factors have no unit!
   if (!evaluator->Evaluate(_sy,sy)) return false;
   if (!evaluator->Evaluate(_sz,sz)) return false;
   if (!evaluator->Evaluate(_rx,rx,aunit)) return false;
   if (!evaluator->Evaluate(_ry,ry,aunit)) return false;
   if (!evaluator->Evaluate(_rz,rz,aunit)) return false;
   if (!evaluator->Evaluate(_dx,dx,lunit)) return false;
   if (!evaluator->Evaluate(_dy,dy,lunit)) return false;
   if (!evaluator->Evaluate(_dz,dz,lunit)) return false;

   G4VSolid* solidPtr = getSolid(nameProcess(solid));

   G4RotationMatrix rot;
   
   rot.rotateX(_rx);
   rot.rotateY(_ry);
   rot.rotateZ(_rz);

   G4ThreeVector trans(_dx,_dy,_dz);

   G4Scale3D scale(_sx,_sy,_sz);

   G4Transform3D transform(rot,trans);
   transform = transform*scale;
          
   new G4ReflectedSolid(nameProcess(name),solidPtr,transform);

   return true;
}

G4bool G4GDMLSolids::rotationRead(const xercesc::DOMElement* const element,G4ThreeVector& vect) {

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

   G4double _x,_y,_z;

   if (!evaluator->Evaluate(_x,x,unit)) return false;
   if (!evaluator->Evaluate(_y,y,unit)) return false;
   if (!evaluator->Evaluate(_z,z,unit)) return false;

   vect.set(_x,_y,_z);

   return true;
}

G4bool G4GDMLSolids::sectionRead(const xercesc::DOMElement* const element,G4ExtrudedSolid::ZSection& section,const G4String& lunit) {

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

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="zPosition"    ) { zPosition     = attribute_value; } else
      if (attribute_name=="xOffset"      ) { xOffset       = attribute_value; } else
      if (attribute_name=="yOffset"      ) { yOffset       = attribute_value; } else
      if (attribute_name=="scalingFactor") { scalingFactor = attribute_value; }
   }

   G4double _zPosition;
   G4double _xOffset;
   G4double _yOffset;
   G4double _scalingFactor;

   if (!evaluator->Evaluate(_zPosition    ,zPosition ,lunit)) return false;
   if (!evaluator->Evaluate(_xOffset      ,xOffset   ,lunit)) return false;
   if (!evaluator->Evaluate(_yOffset      ,yOffset   ,lunit)) return false;
   if (!evaluator->Evaluate(_scalingFactor,scalingFactor)) return false; // scaling factor has no unit!

   section.fZ = _zPosition;
   section.fOffset.set(_xOffset,_yOffset);
   section.fScale = _scalingFactor;

   return true;
}

G4bool G4GDMLSolids::sphereRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit;
   G4String aunit;
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

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name"      ) { name       = attribute_value; } else
      if (attribute_name=="lunit"     ) { lunit      = attribute_value; } else
      if (attribute_name=="aunit"     ) { aunit      = attribute_value; } else
      if (attribute_name=="rmin"      ) { rmin       = attribute_value; } else
      if (attribute_name=="rmax"      ) { rmax       = attribute_value; } else
      if (attribute_name=="startphi"  ) { startphi   = attribute_value; } else
      if (attribute_name=="deltaphi"  ) { deltaphi   = attribute_value; } else
      if (attribute_name=="starttheta") { starttheta = attribute_value; } else
      if (attribute_name=="deltatheta") { deltatheta = attribute_value; }
   }

   G4double _rmin;
   G4double _rmax;
   G4double _startphi;
   G4double _deltaphi;
   G4double _starttheta;
   G4double _deltatheta;

   if (!evaluator->Evaluate(_rmin      ,rmin      ,lunit)) return false;
   if (!evaluator->Evaluate(_rmax      ,rmax      ,lunit)) return false;
   if (!evaluator->Evaluate(_startphi  ,startphi  ,aunit)) return false;
   if (!evaluator->Evaluate(_deltaphi  ,deltaphi  ,aunit)) return false;
   if (!evaluator->Evaluate(_starttheta,starttheta,aunit)) return false;
   if (!evaluator->Evaluate(_deltatheta,deltatheta,aunit)) return false;

   new G4Sphere(nameProcess(name),_rmin,_rmax,_startphi,_deltaphi,_starttheta,_deltatheta);

   return true;
}

G4bool G4GDMLSolids::tessellatedRead(const xercesc::DOMElement* const element) {

   G4String name;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") { name = attribute_value; }
   }
   
   G4TessellatedSolid *tessellated = new G4TessellatedSolid(nameProcess(name));

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="triangular"  ) { if (!triangularRead  (child,tessellated)) return false; } else
      if (tag=="quadrangular") { if (!quadrangularRead(child,tessellated)) return false; }
   }

   tessellated->SetSolidClosed(true);

   return true;
}

G4bool G4GDMLSolids::tetRead(const xercesc::DOMElement* const element) {

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

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name   ") { name    = attribute_value; } else
      if (attribute_name=="vertex1") { vertex1 = attribute_value; } else
      if (attribute_name=="vertex2") { vertex2 = attribute_value; } else
      if (attribute_name=="vertex3") { vertex3 = attribute_value; } else
      if (attribute_name=="vertex4") { vertex4 = attribute_value; }
   }
   
   const G4ThreeVector* ptr1 = define.getPosition(nameProcess(vertex1));
   const G4ThreeVector* ptr2 = define.getPosition(nameProcess(vertex2));
   const G4ThreeVector* ptr3 = define.getPosition(nameProcess(vertex3));
   const G4ThreeVector* ptr4 = define.getPosition(nameProcess(vertex4));

   if (!ptr1 || !ptr2 || !ptr3 || !ptr4) return false;
   
   new G4Tet(nameProcess(name),*ptr1,*ptr2,*ptr3,*ptr4);
   
   return true;
}


G4bool G4GDMLSolids::torusRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit;
   G4String aunit;
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

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name"      ) { name       = attribute_value; } else
      if (attribute_name=="lunit"     ) { lunit      = attribute_value; } else
      if (attribute_name=="aunit"     ) { aunit      = attribute_value; } else
      if (attribute_name=="rmin"      ) { rmin       = attribute_value; } else
      if (attribute_name=="rmax"      ) { rmax       = attribute_value; } else
      if (attribute_name=="rtor"      ) { rtor       = attribute_value; } else
      if (attribute_name=="startphi"  ) { startphi   = attribute_value; } else
      if (attribute_name=="deltaphi"  ) { deltaphi   = attribute_value; }
   }

   G4double _rmin;
   G4double _rmax;
   G4double _rtor;
   G4double _startphi;
   G4double _deltaphi;

   if (!evaluator->Evaluate(_rmin    ,rmin    ,lunit)) return false;
   if (!evaluator->Evaluate(_rmax    ,rmax    ,lunit)) return false;
   if (!evaluator->Evaluate(_rtor    ,rtor    ,lunit)) return false;
   if (!evaluator->Evaluate(_startphi,startphi,aunit)) return false;
   if (!evaluator->Evaluate(_deltaphi,deltaphi,aunit)) return false;

   new G4Torus(nameProcess(name),_rmin,_rmax,_rtor,_startphi,_deltaphi);

   return true;
}

G4bool G4GDMLSolids::trapRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit;
   G4String aunit;
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

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name"  ) { name   = attribute_value; } else
      if (attribute_name=="lunit" ) { lunit  = attribute_value; } else
      if (attribute_name=="aunit" ) { aunit  = attribute_value; } else
      if (attribute_name=="z"     ) { z      = attribute_value; } else
      if (attribute_name=="theta" ) { theta  = attribute_value; } else
      if (attribute_name=="phi"   ) { phi    = attribute_value; } else
      if (attribute_name=="y1"    ) { y1     = attribute_value; } else
      if (attribute_name=="x1"    ) { x1     = attribute_value; } else
      if (attribute_name=="x2"    ) { x2     = attribute_value; } else
      if (attribute_name=="alpha1") { alpha1 = attribute_value; } else
      if (attribute_name=="y2"    ) { y2     = attribute_value; } else
      if (attribute_name=="x3"    ) { x3     = attribute_value; } else
      if (attribute_name=="x4"    ) { x4     = attribute_value; } else
      if (attribute_name=="alpha2") { alpha2 = attribute_value; }
   }

   G4double _z;
   G4double _theta;
   G4double _phi;
   G4double _y1;
   G4double _x1;
   G4double _x2;
   G4double _alpha1;
   G4double _y2;
   G4double _x3;
   G4double _x4;
   G4double _alpha2;

   if (!evaluator->Evaluate(_z     ,z     ,lunit)) return false;
   if (!evaluator->Evaluate(_theta ,theta ,aunit)) return false;
   if (!evaluator->Evaluate(_phi   ,phi   ,aunit)) return false;
   if (!evaluator->Evaluate(_y1    ,y1    ,lunit)) return false;
   if (!evaluator->Evaluate(_x1    ,x1    ,lunit)) return false;
   if (!evaluator->Evaluate(_x2    ,x2    ,lunit)) return false;
   if (!evaluator->Evaluate(_alpha1,alpha1,aunit)) return false;
   if (!evaluator->Evaluate(_y2    ,y2    ,lunit)) return false;
   if (!evaluator->Evaluate(_x3    ,x3    ,lunit)) return false;
   if (!evaluator->Evaluate(_x4    ,x4    ,lunit)) return false;
   if (!evaluator->Evaluate(_alpha2,alpha2,aunit)) return false;

   _z  *= 0.5;
   _y1 *= 0.5;
   _x1 *= 0.5;
   _x2 *= 0.5;
   _y2 *= 0.5;
   _x3 *= 0.5;
   _x4 *= 0.5;

   new G4Trap(nameProcess(name),_z,_theta,_phi,_y1,_x1,_x2,_alpha1,_y2,_x3,_x4,_alpha2);

   return true;
}

G4bool G4GDMLSolids::trdRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit;
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

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name"  ) { name   = attribute_value; } else
      if (attribute_name=="lunit" ) { lunit  = attribute_value; } else
      if (attribute_name=="x1"    ) { x1     = attribute_value; } else
      if (attribute_name=="x2"    ) { x2     = attribute_value; } else
      if (attribute_name=="y1"    ) { y1     = attribute_value; } else
      if (attribute_name=="y2"    ) { y2     = attribute_value; } else
      if (attribute_name=="z"     ) { z      = attribute_value; }
   }

   G4double _x1;
   G4double _x2;
   G4double _y1;
   G4double _y2;
   G4double _z;

   if (!evaluator->Evaluate(_x1,x1,lunit)) return false;
   if (!evaluator->Evaluate(_x2,x2,lunit)) return false;
   if (!evaluator->Evaluate(_y1,y1,lunit)) return false;
   if (!evaluator->Evaluate(_y2,y2,lunit)) return false;
   if (!evaluator->Evaluate(_z ,z ,lunit)) return false;

   _x1 *= 0.5;
   _x2 *= 0.5;
   _y1 *= 0.5;
   _y2 *= 0.5;
   _z  *= 0.5;

   new G4Trd(nameProcess(name),_x1,_x2,_y1,_y2,_z);

   return true;
}

G4bool G4GDMLSolids::triangularRead(const xercesc::DOMElement* const element,G4TessellatedSolid* tessellated) {

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

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="v1"  ) { v1   = attribute_value; } else
      if (attribute_name=="v2"  ) { v2   = attribute_value; } else
      if (attribute_name=="v3"  ) { v3   = attribute_value; } else
      if (attribute_name=="type") { type = attribute_value; }
   }

   const G4ThreeVector* ptr1 = define.getPosition(nameProcess(v1));
   const G4ThreeVector* ptr2 = define.getPosition(nameProcess(v2));
   const G4ThreeVector* ptr3 = define.getPosition(nameProcess(v3));

   if (!ptr1 || !ptr2 || !ptr3) return false;

   tessellated->AddFacet(new G4TriangularFacet(*ptr1,*ptr2,*ptr3,(type=="RELATIVE")?(RELATIVE):(ABSOLUTE)));

   return true;
}

G4bool G4GDMLSolids::tubeRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit;
   G4String aunit;
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

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name"    ) { name     = attribute_value; } else
      if (attribute_name=="lunit"   ) { lunit    = attribute_value; } else
      if (attribute_name=="aunit"   ) { aunit    = attribute_value; } else
      if (attribute_name=="rmin"    ) { rmin     = attribute_value; } else
      if (attribute_name=="rmax"    ) { rmax     = attribute_value; } else
      if (attribute_name=="z"       ) { z        = attribute_value; } else
      if (attribute_name=="startphi") { startphi = attribute_value; } else
      if (attribute_name=="deltaphi") { deltaphi = attribute_value; }
   }

   G4double _rmin;
   G4double _rmax;
   G4double _z;
   G4double _startphi;
   G4double _deltaphi;

   if (!evaluator->Evaluate(_rmin    ,rmin    ,lunit)) return false;
   if (!evaluator->Evaluate(_rmax    ,rmax    ,lunit)) return false;
   if (!evaluator->Evaluate(_z       ,z       ,lunit)) return false;
   if (!evaluator->Evaluate(_startphi,startphi,aunit)) return false;
   if (!evaluator->Evaluate(_deltaphi,deltaphi,aunit)) return false;

   _z *= 0.5;

   new G4Tubs(nameProcess(name),_rmin,_rmax,_z,_startphi,_deltaphi);

   return true;
}

G4bool G4GDMLSolids::twoDimVertexRead(const xercesc::DOMElement* const element,G4TwoVector& vec2D,const G4String& lunit) {

   G4String x;
   G4String y;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="x") { x = attribute_value; } else
      if (attribute_name=="y") { y = attribute_value; }
   }

   G4double _x;
   G4double _y;

   if (!evaluator->Evaluate(_x,x,lunit)) return false;
   if (!evaluator->Evaluate(_y,y,lunit)) return false;

   vec2D.set(_x,_y);

   return true;
}

G4bool G4GDMLSolids::xtruRead(const xercesc::DOMElement* const element) {

   G4String name;
   G4String lunit;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name" ) { name  = attribute_value; } else
      if (attribute_name=="lunit") { lunit = attribute_value; }
   }

   std::vector<G4TwoVector> twoDimVertexList;
   std::vector<G4ExtrudedSolid::ZSection> sectionList;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="twoDimVertex") {
      
         G4TwoVector temp;
	 if (!twoDimVertexRead(child,temp,lunit)) return false;
         twoDimVertexList.push_back(temp);
      } else
      if (tag=="section") {

         G4ExtrudedSolid::ZSection temp(0,G4TwoVector(),0);
	 if (!sectionRead(child,temp,lunit)) return false;
	 sectionList.push_back(temp);      
      }
   }

   new G4ExtrudedSolid(nameProcess(name),twoDimVertexList,sectionList);

   return true;
}

G4bool G4GDMLSolids::zplaneRead(const xercesc::DOMElement* const element,zplaneType& zplane,const G4String& lunit) {

   G4String rmin;
   G4String rmax;
   G4String z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const G4String attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const G4String attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="rmin" ) { rmin  = attribute_value; } else
      if (attribute_name=="rmax" ) { rmax  = attribute_value; } else
      if (attribute_name=="z"    ) { z     = attribute_value; }
   }

   G4double _rmin;
   G4double _rmax;
   G4double _z;

   if (!evaluator->Evaluate(_rmin,rmin,lunit)) return false;
   if (!evaluator->Evaluate(_rmax,rmax,lunit)) return false;
   if (!evaluator->Evaluate(_z   ,z   ,lunit)) return false;

   zplane.rmin = _rmin;
   zplane.rmax = _rmax;
   zplane.z    = _z;

   return true;
}

G4bool G4GDMLSolids::Read(const xercesc::DOMElement* const element,G4GDMLEvaluator *evalPtr,const G4String& file0) {

   evaluator = evalPtr;
   
   file = file0;

   for (xercesc::DOMNode* iter = element->getFirstChild();iter != 0;iter = iter->getNextSibling()) {

      if (iter->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;

      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(iter);

      const G4String tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="box"           ) { if (!boxRead           (child)) return false; } else
      if (tag=="cone"          ) { if (!coneRead          (child)) return false; } else
      if (tag=="ellipsoid"     ) { if (!ellipsoidRead     (child)) return false; } else
      if (tag=="eltube"        ) { if (!eltubeRead        (child)) return false; } else
      if (tag=="hype"          ) { if (!hypeRead          (child)) return false; } else
      if (tag=="loop"          ) { if (!loopRead          (child)) return false; } else
      if (tag=="orb"           ) { if (!orbRead           (child)) return false; } else
      if (tag=="para"          ) { if (!paraRead          (child)) return false; } else
      if (tag=="polycone"      ) { if (!polyconeRead      (child)) return false; } else
      if (tag=="polyhedra"     ) { if (!polyhedraRead     (child)) return false; } else
      if (tag=="sphere"        ) { if (!sphereRead        (child)) return false; } else
      if (tag=="reflectedSolid") { if (!reflectedSolidRead(child)) return false; } else
      if (tag=="tessellated"   ) { if (!tessellatedRead   (child)) return false; } else
      if (tag=="tet"           ) { if (!tetRead           (child)) return false; } else
      if (tag=="torus"         ) { if (!torusRead         (child)) return false; } else
      if (tag=="trap"          ) { if (!trapRead          (child)) return false; } else
      if (tag=="trd"           ) { if (!trdRead           (child)) return false; } else
      if (tag=="tube"          ) { if (!tubeRead          (child)) return false; } else
      if (tag=="xtru"          ) { if (!xtruRead          (child)) return false; } else
      if (tag=="intersection"  ) { if (!booleanRead(child,INTERSECTION)) return false; } else // Parse union,intersection and subtraction with the same function!
      if (tag=="subtraction"   ) { if (!booleanRead(child,SUBTRACTION )) return false; } else
      if (tag=="union"         ) { if (!booleanRead(child,UNION       )) return false; } else
      G4Exception("GDML: Unknown tag in solids: "+tag);
   }

   return true;
}

G4VSolid* G4GDMLSolids::getSolid(const G4String& ref) const {

   G4VSolid *solidPtr = G4SolidStore::GetInstance()->GetSolid(ref,false);

   if (!solidPtr) G4Exception("GDML: Referenced solid '"+ref+"' not found!");

   return solidPtr;
}
