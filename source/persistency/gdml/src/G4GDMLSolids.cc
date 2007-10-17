#include "G4GDMLSolids.hh"

G4GDMLSolids::G4GDMLSolids() {

   evaluator = G4GDMLEvaluator::GetInstance();
}

bool G4GDMLSolids::booleanRead(const xercesc::DOMElement* const element,const BooleanOp op) {

   std::string name;
   std::string first_ref;
   std::string second_ref;

   G4ThreeVector position;
   G4ThreeVector rotation;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="name") { name = attribute_value; } else
      {
      }
   }

   const xercesc::DOMNodeList* const children = element->getChildNodes();
   XMLSize_t elementCount = children->getLength();

   for (XMLSize_t element_index=0;element_index<elementCount;element_index++) {

      xercesc::DOMNode* element_node = children->item(element_index);
      
      if (element_node->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;
   
      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(element_node);   

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="first"   ) { if (!refRead     (child,first_ref )) return false; } else
      if (tag=="second"  ) { if (!refRead     (child,second_ref)) return false; } else
      if (tag=="position") { if (!positionRead(child,position  )) return false; } else
      if (tag=="rotation") { if (!rotationRead(child,rotation  )) return false; } else
      {
         std::cout << std::endl;
         std::cout << "GDML ERROR! Unsupported tag in boolean solid '" << name << "': " << tag << std::endl;
         std::cout << std::endl;
         return false;
      }
   }

   G4VSolid *FirstSolid = G4SolidStore::GetInstance()->GetSolid(first_ref);

   if (!FirstSolid) {

      std::cout << std::endl;   
      std::cout << "GDML ERROR! Referenced solid '" << first_ref << "' in boolean solid '" << name << "' was not found!" << std::endl;   
      std::cout << std::endl;   
      return false;
   }

   G4VSolid *SecondSolid = G4SolidStore::GetInstance()->GetSolid(second_ref);

   if (!SecondSolid) {

      std::cout << std::endl;   
      std::cout << "GDML ERROR! Referenced solid '" << second_ref << "' in boolean solid '" << name << "' was not found!" << std::endl;   
      std::cout << std::endl;   
      return false;
   }

   G4RotationMatrix rot;

   rot.rotateX(rotation.x());
   rot.rotateY(rotation.y());
   rot.rotateZ(rotation.z());

   G4Transform3D transform(rot,position);

   if (op==UNION       ) { new G4UnionSolid       (name,FirstSolid,SecondSolid,transform); } else
   if (op==SUBTRACTION ) { new G4SubtractionSolid (name,FirstSolid,SecondSolid,transform); } else
   if (op==INTERSECTION) { new G4IntersectionSolid(name,FirstSolid,SecondSolid,transform); }

   return true;
}

bool G4GDMLSolids::boxRead(const xercesc::DOMElement* const element) {

   std::string lunit;
   std::string aunit;
   std::string name;
   std::string x;
   std::string y;
   std::string z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="lunit") { lunit = attribute_value; } else
      if (attribute_name=="aunit") { aunit = attribute_value; } else
      if (attribute_name=="name" ) { name  = attribute_value; } else
      if (attribute_name=="x"    ) { x     = attribute_value; } else
      if (attribute_name=="y"    ) { y     = attribute_value; } else
      if (attribute_name=="z"    ) { z     = attribute_value; } else
      {
      }
   }

   double _x;
   double _y;
   double _z;

   if (!evaluator->Evaluate(_x,x,lunit)) return false;
   if (!evaluator->Evaluate(_y,y,lunit)) return false;
   if (!evaluator->Evaluate(_z,z,lunit)) return false;

   _x *= 0.5;
   _y *= 0.5;
   _z *= 0.5;

   new G4Box(name,_x,_y,_z);

   return true;
}

bool G4GDMLSolids::coneRead(const xercesc::DOMElement* const element) {

   std::string lunit;
   std::string aunit;
   std::string name;
   std::string rmin1;
   std::string rmax1;
   std::string rmin2;
   std::string rmax2;
   std::string z;
   std::string startphi;
   std::string deltaphi;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="lunit"   ) { lunit    = attribute_value; } else
      if (attribute_name=="aunit"   ) { aunit    = attribute_value; } else
      if (attribute_name=="name"    ) { name     = attribute_value; } else
      if (attribute_name=="rmin1"   ) { rmin1    = attribute_value; } else
      if (attribute_name=="rmax1"   ) { rmax1    = attribute_value; } else
      if (attribute_name=="rmin2"   ) { rmin2    = attribute_value; } else
      if (attribute_name=="rmax2"   ) { rmax2    = attribute_value; } else
      if (attribute_name=="z"       ) { z        = attribute_value; } else
      if (attribute_name=="startphi") { startphi = attribute_value; } else
      if (attribute_name=="deltaphi") { deltaphi = attribute_value; } else
      {
      }
   }

   double _rmin1;
   double _rmax1;
   double _rmin2;
   double _rmax2;
   double _z;
   double _startphi;
   double _deltaphi;

   if (!evaluator->Evaluate(_rmin1   ,rmin1   ,lunit)) return false;
   if (!evaluator->Evaluate(_rmax1   ,rmax1   ,lunit)) return false;
   if (!evaluator->Evaluate(_rmin2   ,rmin2   ,lunit)) return false;
   if (!evaluator->Evaluate(_rmax2   ,rmax2   ,lunit)) return false;
   if (!evaluator->Evaluate(_z       ,z       ,lunit)) return false;
   if (!evaluator->Evaluate(_startphi,startphi,aunit)) return false;
   if (!evaluator->Evaluate(_deltaphi,deltaphi,aunit)) return false;

   _z *= 0.5;

   new G4Cons(name,_rmin1,_rmax1,_rmin2,_rmax2,_z,_startphi,_deltaphi);

   return true;
}

bool G4GDMLSolids::ellipsoidRead(const xercesc::DOMElement* const element) {

   std::string lunit;
   std::string aunit;
   std::string name;
   std::string ax;
   std::string by;
   std::string cz;
   std::string zcut1; 
   std::string zcut2; 

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="lunit") { lunit = attribute_value; } else
      if (attribute_name=="aunit") { aunit = attribute_value; } else
      if (attribute_name=="name" ) { name  = attribute_value; } else
      if (attribute_name=="ax"   ) { ax    = attribute_value; } else
      if (attribute_name=="by"   ) { by    = attribute_value; } else
      if (attribute_name=="cz"   ) { cz    = attribute_value; } else
      if (attribute_name=="zcut1") { zcut1 = attribute_value; } else
      if (attribute_name=="zcut2") { zcut2 = attribute_value; } else
      {
      }
   }

   double _ax;
   double _by;
   double _cz;
   double _zcut1; 
   double _zcut2; 

   if (!evaluator->Evaluate(_ax   ,ax   ,lunit)) return false;
   if (!evaluator->Evaluate(_by   ,by   ,lunit)) return false;
   if (!evaluator->Evaluate(_cz   ,cz   ,lunit)) return false;
   if (!evaluator->Evaluate(_zcut1,zcut1,lunit)) return false;
   if (!evaluator->Evaluate(_zcut2,zcut2,lunit)) return false;

   new G4Ellipsoid(name,_ax,_by,_cz,_zcut1,_zcut2);

   return true;
}

bool G4GDMLSolids::hypeRead(const xercesc::DOMElement* const element) {

   std::string lunit;
   std::string aunit;
   std::string name;
   std::string rmin;
   std::string rmax;
   std::string inst;
   std::string outst;
   std::string z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="lunit") { lunit = attribute_value; } else
      if (attribute_name=="aunit") { aunit = attribute_value; } else
      if (attribute_name=="name" ) { name  = attribute_value; } else
      if (attribute_name=="rmin" ) { rmin  = attribute_value; } else
      if (attribute_name=="rmax" ) { rmax  = attribute_value; } else
      if (attribute_name=="inst" ) { inst  = attribute_value; } else
      if (attribute_name=="outst") { outst = attribute_value; } else
      if (attribute_name=="z"    ) { z     = attribute_value; } else
      {
      }
   }

   double _rmin;
   double _rmax;
   double _inst;
   double _outst;
   double _z;

   if (!evaluator->Evaluate(_rmin ,rmin ,lunit)) return false;
   if (!evaluator->Evaluate(_rmax ,rmax ,lunit)) return false;
   if (!evaluator->Evaluate(_inst ,inst ,aunit)) return false;
   if (!evaluator->Evaluate(_outst,outst,aunit)) return false;
   if (!evaluator->Evaluate(_z    ,z    ,lunit)) return false;

   _z *= 0.5;
   
   new G4Hype(name,_rmin,_rmax,_inst,_outst,_z);

   return true;
}

bool G4GDMLSolids::orbRead(const xercesc::DOMElement* const element) {

   std::string lunit;
   std::string aunit;
   std::string name;
   std::string r;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="lunit") { lunit = attribute_value; } else
      if (attribute_name=="aunit") { aunit = attribute_value; } else
      if (attribute_name=="name" ) { name  = attribute_value; } else
      if (attribute_name=="r"    ) { r     = attribute_value; } else
      {
      }
   }

   double _r;

   if (!evaluator->Evaluate(_r,r,lunit)) return false;

   new G4Orb(name,_r);

   return true;
}

bool G4GDMLSolids::paraRead(const xercesc::DOMElement* const element) {

   std::string lunit;
   std::string aunit;
   std::string name;
   std::string x;
   std::string y;
   std::string z;
   std::string alpha;
   std::string theta;
   std::string phi;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="lunit") { lunit = attribute_value; } else
      if (attribute_name=="aunit") { aunit = attribute_value; } else
      if (attribute_name=="name" ) { name  = attribute_value; } else
      if (attribute_name=="x"    ) { x     = attribute_value; } else
      if (attribute_name=="y"    ) { y     = attribute_value; } else
      if (attribute_name=="z"    ) { z     = attribute_value; } else
      if (attribute_name=="alpha") { alpha = attribute_value; } else
      if (attribute_name=="theta") { theta = attribute_value; } else
      if (attribute_name=="phi"  ) { phi   = attribute_value; } else
      {
      }
   }

   double _x;
   double _y;
   double _z;
   double _alpha;
   double _theta;
   double _phi;

   if (!evaluator->Evaluate(_x    ,x    ,lunit)) return false;
   if (!evaluator->Evaluate(_y    ,y    ,lunit)) return false;
   if (!evaluator->Evaluate(_z    ,z    ,lunit)) return false;
   if (!evaluator->Evaluate(_alpha,alpha,aunit)) return false;
   if (!evaluator->Evaluate(_theta,theta,aunit)) return false;
   if (!evaluator->Evaluate(_phi  ,phi  ,aunit)) return false;

   _x *= 0.5;
   _y *= 0.5;
   _z *= 0.5;

   new G4Para(name,_x,_y,_z,_alpha,_theta,_phi);

   return true;
}

bool G4GDMLSolids::polyconeRead(const xercesc::DOMElement* const element) {

   std::string lunit;
   std::string aunit;
   std::string name;
   std::string startphi;
   std::string deltaphi;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="lunit"   ) { lunit    = attribute_value; } else
      if (attribute_name=="aunit"   ) { aunit    = attribute_value; } else
      if (attribute_name=="name"    ) { name     = attribute_value; } else
      if (attribute_name=="startphi") { startphi = attribute_value; } else
      if (attribute_name=="deltaphi") { deltaphi = attribute_value; } else
      {
      }
   }

   double _startphi;
   double _deltaphi;

   if (!evaluator->Evaluate(_startphi,startphi,aunit)) return false;
   if (!evaluator->Evaluate(_deltaphi,deltaphi,aunit)) return false;

   std::vector<zplaneType> zplaneList;

   const xercesc::DOMNodeList* const children = element->getChildNodes();
   XMLSize_t elementCount = children->getLength();

   for (XMLSize_t element_index=0;element_index<elementCount;element_index++) {

      xercesc::DOMNode* element_node = children->item(element_index);
      
      if (element_node->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;
   
      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(element_node);   

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="zplane") {

         zplaneType zplane;      
      
         if (!zplaneRead(child,zplane,lunit)) return false;

         zplaneList.push_back(zplane);
      }
   }

   int numZPlanes = zplaneList.size();

   double *rmin_array = new double[numZPlanes];
   double *rmax_array = new double[numZPlanes];
   double* z_array = new double[numZPlanes];

   for (int i=0;i<numZPlanes;i++) {
   
      rmin_array[i] = zplaneList[i].rmin;
      rmax_array[i] = zplaneList[i].rmax;
      z_array[i]    = zplaneList[i].z;
   }

   new G4Polycone(name,_startphi,_deltaphi,numZPlanes,z_array,rmin_array,rmax_array);

   return true;
}

bool G4GDMLSolids::polyhedraRead(const xercesc::DOMElement* const element) {

   std::string lunit;
   std::string aunit;
   std::string name;
   std::string startphi;
   std::string totalphi;
   std::string numsides;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="lunit"   ) { lunit    = attribute_value; } else
      if (attribute_name=="aunit"   ) { aunit    = attribute_value; } else
      if (attribute_name=="name"    ) { name     = attribute_value; } else
      if (attribute_name=="startphi") { startphi = attribute_value; } else
      if (attribute_name=="totalphi") { totalphi = attribute_value; } else
      if (attribute_name=="numsides") { numsides = attribute_value; } else
      {
      }
   }

   double _startphi;
   double _totalphi;
   double _numsides;

   if (!evaluator->Evaluate(_startphi,startphi,aunit)) return false;
   if (!evaluator->Evaluate(_totalphi,totalphi,aunit)) return false;

   if (!evaluator->Evaluate(_numsides,numsides)) return false;

   std::vector<zplaneType> zplaneList;

   const xercesc::DOMNodeList* const children = element->getChildNodes();
   XMLSize_t elementCount = children->getLength();

   for (XMLSize_t element_index=0;element_index<elementCount;element_index++) {

      xercesc::DOMNode* element_node = children->item(element_index);
      
      if (element_node->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;
   
      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(element_node);   

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="zplane") {

         zplaneType zplane;      
      
         if (!zplaneRead(child,zplane,lunit)) return false;

         zplaneList.push_back(zplane);
      }
   }

   int numZPlanes = zplaneList.size();

   double *rmin_array = new double[numZPlanes];
   double *rmax_array = new double[numZPlanes];
   double* z_array = new double[numZPlanes];

   for (int i=0;i<numZPlanes;i++) {
   
      rmin_array[i] = zplaneList[i].rmin;
      rmax_array[i] = zplaneList[i].rmax;
      z_array[i]    = zplaneList[i].z;
   }

   new G4Polyhedra(name,_startphi,_totalphi,(int)_numsides,numZPlanes,z_array,rmin_array,rmax_array);

   return true;
}

bool G4GDMLSolids::positionRead(const xercesc::DOMElement* const element,G4ThreeVector& vect) {

   std::string unit;
   std::string name;
   std::string x;
   std::string y;
   std::string z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="unit") { unit = attribute_value; } else
      if (attribute_name=="name") { name = attribute_value; } else
      if (attribute_name=="x"   ) { x    = attribute_value; } else
      if (attribute_name=="y"   ) { y    = attribute_value; } else
      if (attribute_name=="z"   ) { z    = attribute_value; } else
      {
      }
   }

   double _x;
   double _y;
   double _z;

   if (!evaluator->Evaluate(_x,x,unit)) return false;
   if (!evaluator->Evaluate(_y,y,unit)) return false;
   if (!evaluator->Evaluate(_z,z,unit)) return false;

   vect.set(_x,_y,_z);

   return true;
}

bool G4GDMLSolids::quadrangularRead(const xercesc::DOMElement* const element) {

   return true;
}

bool G4GDMLSolids::refRead(const xercesc::DOMElement* const element,std::string& ref) {

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   if (attributeCount != 1) return false;

   xercesc::DOMNode* attribute_node = attributes->item(0);

   if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) return true;

   const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

   const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
   const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

   if (attribute_name=="ref") { ref = attribute_value; } else
   {
   }

   return true;
}

bool G4GDMLSolids::rotationRead(const xercesc::DOMElement* const element,G4ThreeVector& vect) {

   std::string unit;
   std::string name;
   std::string x;
   std::string y;
   std::string z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="unit") { unit = attribute_value; } else
      if (attribute_name=="name") { name = attribute_value; } else
      if (attribute_name=="x"   ) { x    = attribute_value; } else
      if (attribute_name=="y"   ) { y    = attribute_value; } else
      if (attribute_name=="z"   ) { z    = attribute_value; } else
      {
      }
   }

   double _x;
   double _y;
   double _z;

   if (!evaluator->Evaluate(_x,x,unit)) return false;
   if (!evaluator->Evaluate(_y,y,unit)) return false;
   if (!evaluator->Evaluate(_z,z,unit)) return false;

   vect.set(_x,_y,_z);

   return true;
}

bool G4GDMLSolids::sectionRead(const xercesc::DOMElement* const element,G4ExtrudedSolid::ZSection& section,const std::string& lunit) {

   std::string zPosition;
   std::string xOffset;
   std::string yOffset;
   std::string scalingFactor;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="zPosition"    ) { zPosition     = attribute_value; } else
      if (attribute_name=="xOffset"      ) { xOffset       = attribute_value; } else
      if (attribute_name=="yOffset"      ) { yOffset       = attribute_value; } else
      if (attribute_name=="scalingFactor") { scalingFactor = attribute_value; } else
      {
      }
   }

   double _zPosition;
   double _xOffset;
   double _yOffset;
   double _scalingFactor;

   if (!evaluator->Evaluate(_zPosition    ,zPosition    ,lunit)) return false;
   if (!evaluator->Evaluate(_xOffset      ,xOffset      ,lunit)) return false;
   if (!evaluator->Evaluate(_yOffset      ,yOffset      ,lunit)) return false;
   if (!evaluator->Evaluate(_scalingFactor,scalingFactor,lunit)) return false;

   section.fZ = _zPosition;
   section.fOffset.set(_xOffset,_yOffset);
   section.fScale = _scalingFactor;

   return true;
}

bool G4GDMLSolids::sphereRead(const xercesc::DOMElement* const element) {

   std::string lunit;
   std::string aunit;
   std::string name;
   std::string rmin;
   std::string rmax;
   std::string startphi;
   std::string deltaphi;
   std::string starttheta;
   std::string deltatheta;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="lunit"     ) { lunit      = attribute_value; } else
      if (attribute_name=="aunit"     ) { aunit      = attribute_value; } else
      if (attribute_name=="name"      ) { name       = attribute_value; } else
      if (attribute_name=="rmin"      ) { rmin       = attribute_value; } else
      if (attribute_name=="rmax"      ) { rmax       = attribute_value; } else
      if (attribute_name=="startphi"  ) { startphi   = attribute_value; } else
      if (attribute_name=="deltaphi"  ) { deltaphi   = attribute_value; } else
      if (attribute_name=="starttheta") { starttheta = attribute_value; } else
      if (attribute_name=="deltatheta") { deltatheta = attribute_value; } else
      {
      }
   }

   double _rmin;
   double _rmax;
   double _startphi;
   double _deltaphi;
   double _starttheta;
   double _deltatheta;

   if (!evaluator->Evaluate(_rmin      ,rmin      ,lunit)) return false;
   if (!evaluator->Evaluate(_rmax      ,rmax      ,lunit)) return false;
   if (!evaluator->Evaluate(_startphi  ,startphi  ,aunit)) return false;
   if (!evaluator->Evaluate(_deltaphi  ,deltaphi  ,aunit)) return false;
   if (!evaluator->Evaluate(_starttheta,starttheta,aunit)) return false;
   if (!evaluator->Evaluate(_deltatheta,deltatheta,aunit)) return false;

   new G4Sphere(name,_rmin,_rmax,_startphi,_deltaphi,_starttheta,_deltatheta);

   return true;
}

bool G4GDMLSolids::tessellatedRead(const xercesc::DOMElement* const element) {

   std::string lunit;
   std::string aunit;
   std::string name;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="lunit"     ) { lunit      = attribute_value; } else
      if (attribute_name=="aunit"     ) { aunit      = attribute_value; } else
      if (attribute_name=="name"      ) { name       = attribute_value; } else
      {
      }
   }
   
   G4TessellatedSolid *tessellated = new G4TessellatedSolid(name);

   return true;
}

bool G4GDMLSolids::torusRead(const xercesc::DOMElement* const element) {

   std::string lunit;
   std::string aunit;
   std::string name;
   std::string rmin;
   std::string rmax;
   std::string rtor;
   std::string startphi;
   std::string deltaphi;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="lunit"     ) { lunit      = attribute_value; } else
      if (attribute_name=="aunit"     ) { aunit      = attribute_value; } else
      if (attribute_name=="name"      ) { name       = attribute_value; } else
      if (attribute_name=="rmin"      ) { rmin       = attribute_value; } else
      if (attribute_name=="rmax"      ) { rmax       = attribute_value; } else
      if (attribute_name=="rtor"      ) { rtor       = attribute_value; } else
      if (attribute_name=="startphi"  ) { startphi   = attribute_value; } else
      if (attribute_name=="deltaphi"  ) { deltaphi   = attribute_value; } else
      {
      }
   }

   double _rmin;
   double _rmax;
   double _rtor;
   double _startphi;
   double _deltaphi;

   if (!evaluator->Evaluate(_rmin    ,rmin    ,lunit)) return false;
   if (!evaluator->Evaluate(_rmax    ,rmax    ,lunit)) return false;
   if (!evaluator->Evaluate(_rtor    ,rtor    ,lunit)) return false;
   if (!evaluator->Evaluate(_startphi,startphi,aunit)) return false;
   if (!evaluator->Evaluate(_deltaphi,deltaphi,aunit)) return false;

   new G4Torus(name,_rmin,_rmax,_rtor,_startphi,_deltaphi);

   return true;
}

bool G4GDMLSolids::trapRead(const xercesc::DOMElement* const element) {

   std::string lunit;
   std::string aunit;
   std::string name;
   std::string z;
   std::string theta;
   std::string phi;
   std::string y1;
   std::string x1;
   std::string x2;
   std::string alpha1;
   std::string y2;
   std::string x3;
   std::string x4;
   std::string alpha2;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="lunit" ) { lunit  = attribute_value; } else
      if (attribute_name=="aunit" ) { aunit  = attribute_value; } else
      if (attribute_name=="name"  ) { name   = attribute_value; } else
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
      if (attribute_name=="alpha2") { alpha2 = attribute_value; } else
      {
      }
   }

   double _z;
   double _theta;
   double _phi;
   double _y1;
   double _x1;
   double _x2;
   double _alpha1;
   double _y2;
   double _x3;
   double _x4;
   double _alpha2;

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

   new G4Trap(name,_z,_theta,_phi,_y1,_x1,_x2,_alpha1,_y2,_x3,_x4,_alpha2);

   return true;
}

bool G4GDMLSolids::trdRead(const xercesc::DOMElement* const element) {

   std::string lunit;
   std::string aunit;
   std::string name;
   std::string x1;
   std::string x2;
   std::string y1;
   std::string y2;
   std::string z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="lunit" ) { lunit  = attribute_value; } else
      if (attribute_name=="aunit" ) { aunit  = attribute_value; } else
      if (attribute_name=="name"  ) { name   = attribute_value; } else
      if (attribute_name=="x1"    ) { x1     = attribute_value; } else
      if (attribute_name=="x2"    ) { x2     = attribute_value; } else
      if (attribute_name=="y1"    ) { y1     = attribute_value; } else
      if (attribute_name=="y2"    ) { y2     = attribute_value; } else
      if (attribute_name=="z"     ) { z      = attribute_value; } else
      {
      }
   }

   double _x1;
   double _x2;
   double _y1;
   double _y2;
   double _z;

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

   new G4Trd(name,_x1,_x2,_y1,_y2,_z);

   return true;
}

bool G4GDMLSolids::triangularRead(const xercesc::DOMElement* const element) {

   return true;
}

bool G4GDMLSolids::tubeRead(const xercesc::DOMElement* const element) {

   std::string lunit;
   std::string aunit;
   std::string name;
   std::string rmin;
   std::string rmax;
   std::string z;
   std::string startphi;
   std::string deltaphi;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="lunit"   ) { lunit    = attribute_value; } else
      if (attribute_name=="aunit"   ) { aunit    = attribute_value; } else
      if (attribute_name=="name"    ) { name     = attribute_value; } else
      if (attribute_name=="rmin"    ) { rmin     = attribute_value; } else
      if (attribute_name=="rmax"    ) { rmax     = attribute_value; } else
      if (attribute_name=="z"       ) { z        = attribute_value; } else
      if (attribute_name=="startphi") { startphi = attribute_value; } else
      if (attribute_name=="deltaphi") { deltaphi = attribute_value; } else
      {
      }
   }

   double _rmin;
   double _rmax;
   double _z;
   double _startphi;
   double _deltaphi;

   if (!evaluator->Evaluate(_rmin    ,rmin    ,lunit)) return false;
   if (!evaluator->Evaluate(_rmax    ,rmax    ,lunit)) return false;
   if (!evaluator->Evaluate(_z       ,z       ,lunit)) return false;
   if (!evaluator->Evaluate(_startphi,startphi,aunit)) return false;
   if (!evaluator->Evaluate(_deltaphi,deltaphi,aunit)) return false;

   _z  *= 0.5;

   new G4Tubs(name,_rmin,_rmax,_z,_startphi,_deltaphi);

   return true;
}

bool G4GDMLSolids::twoDimVertexRead(const xercesc::DOMElement* const element,G4TwoVector& vec2D,const std::string& lunit) {

   std::string x;
   std::string y;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="x") { x = attribute_value; } else
      if (attribute_name=="y") { y = attribute_value; } else
      {
      }
   }

   double _x;
   double _y;

   if (!evaluator->Evaluate(_x,x,lunit)) return false;
   if (!evaluator->Evaluate(_y,y,lunit)) return false;

   vec2D.set(_x,_y);

   return true;
}

bool G4GDMLSolids::xtruRead(const xercesc::DOMElement* const element) {

   std::string lunit;
   std::string aunit;
   std::string name;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="lunit"   ) { lunit    = attribute_value; } else
      if (attribute_name=="aunit"   ) { aunit    = attribute_value; } else
      if (attribute_name=="name"    ) { name     = attribute_value; } else
      {
      }
   }

   std::vector<G4TwoVector> twoDimVertexList;
   std::vector<G4ExtrudedSolid::ZSection> sectionList;

   const xercesc::DOMNodeList* const children = element->getChildNodes();
   XMLSize_t elementCount = children->getLength();

   for (XMLSize_t element_index=0;element_index<elementCount;element_index++) {

      xercesc::DOMNode* element_node = children->item(element_index);
      
      if (element_node->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;
   
      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(element_node);   

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

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

   new G4ExtrudedSolid(name,twoDimVertexList,sectionList);

   return true;
}

bool G4GDMLSolids::zplaneRead(const xercesc::DOMElement* const element,zplaneType& zplane,const std::string& lunit) {

   std::string rmin;
   std::string rmax;
   std::string z;

   const xercesc::DOMNamedNodeMap* const attributes = element->getAttributes();
   XMLSize_t attributeCount = attributes->getLength();

   for (XMLSize_t attribute_index=0;attribute_index<attributeCount;attribute_index++) {

      xercesc::DOMNode* attribute_node = attributes->item(attribute_index);

      if (attribute_node->getNodeType() != xercesc::DOMNode::ATTRIBUTE_NODE) continue;

      const xercesc::DOMAttr* const attribute = dynamic_cast<xercesc::DOMAttr*>(attribute_node);   

      const std::string attribute_name  = xercesc::XMLString::transcode(attribute->getName());
      const std::string attribute_value = xercesc::XMLString::transcode(attribute->getValue());

      if (attribute_name=="rmin" ) { rmin  = attribute_value; } else
      if (attribute_name=="rmax" ) { rmax  = attribute_value; } else
      if (attribute_name=="z"    ) { z     = attribute_value; } else
      {
      }
   }

   double _rmin;
   double _rmax;
   double _z;

   if (!evaluator->Evaluate(_rmin,rmin,lunit)) return false;
   if (!evaluator->Evaluate(_rmax,rmax,lunit)) return false;
   if (!evaluator->Evaluate(_z   ,z   ,lunit)) return false;

   zplane.rmin = _rmin;
   zplane.rmax = _rmax;
   zplane.z    = _z;

   return true;
}

bool G4GDMLSolids::Read(const xercesc::DOMElement* const element) {

   const xercesc::DOMNodeList* const children = element->getChildNodes();
   XMLSize_t elementCount = children->getLength();

   for (XMLSize_t element_index=0;element_index<elementCount;element_index++) {

      xercesc::DOMNode* element_node = children->item(element_index);
      
      if (element_node->getNodeType() != xercesc::DOMNode::ELEMENT_NODE) continue;
   
      const xercesc::DOMElement* const child = dynamic_cast<xercesc::DOMElement*>(element_node);   

      const std::string tag = xercesc::XMLString::transcode(child->getTagName());

      if (tag=="box"        ) { if (!boxRead        (child)) return false; } else
      if (tag=="cone"       ) { if (!coneRead       (child)) return false; } else
      if (tag=="ellipsoid"  ) { if (!ellipsoidRead  (child)) return false; } else
      if (tag=="hype"       ) { if (!hypeRead       (child)) return false; } else
      if (tag=="orb"        ) { if (!orbRead        (child)) return false; } else
      if (tag=="para"       ) { if (!paraRead       (child)) return false; } else
      if (tag=="polycone"   ) { if (!polyconeRead   (child)) return false; } else
      if (tag=="polyhedra"  ) { if (!polyhedraRead  (child)) return false; } else
      if (tag=="sphere"     ) { if (!sphereRead     (child)) return false; } else
      if (tag=="tessellated") { if (!tessellatedRead(child)) return false; } else
      if (tag=="torus"      ) { if (!torusRead      (child)) return false; } else
      if (tag=="trap"       ) { if (!trapRead       (child)) return false; } else
      if (tag=="trd"        ) { if (!trdRead        (child)) return false; } else
      if (tag=="tube"       ) { if (!tubeRead       (child)) return false; } else
      if (tag=="xtru"       ) { if (!xtruRead       (child)) return false; } else
      if (tag=="intersection") { if (!booleanRead(child,INTERSECTION)) return false; } else // Parse union,intersection and subtraction with the same function!
      if (tag=="subtraction" ) { if (!booleanRead(child,SUBTRACTION )) return false; } else
      if (tag=="union"       ) { if (!booleanRead(child,UNION       )) return false; } else
      {
         std::cout << std::endl;
	 std::cout << "GDML ERROR! Unsupported tag in solids: " << tag << std::endl;
         std::cout << std::endl;
         return false;
      }
   }

   return true;
}
