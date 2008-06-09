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
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLWriteDefine.hh"

const G4double G4GDMLWriteDefine::kRelativePrecision = DBL_EPSILON;
const G4double G4GDMLWriteDefine::kAngularPrecision = 1.0e-6*degree;
const G4double G4GDMLWriteDefine::kLinearPrecision = 1.0e-15*meter;

G4ThreeVector G4GDMLWriteDefine::getAngles(const G4RotationMatrix& mat) {

   G4double x,y,z;

   G4double cosb = sqrt(mat.xx()*mat.xx()+mat.yx()*mat.yx());

   if (cosb > kRelativePrecision) {

      x = atan2(mat.zy(),mat.zz());
      y = atan2(-mat.zx(),cosb);
      z = atan2(mat.yx(),mat.xx());
   } else {

      x = atan2(-mat.yz(),mat.yy());
      y = atan2(-mat.zx(),cosb);
      z = 0.0;
   }

   return G4ThreeVector(x,y,z);
}

void G4GDMLWriteDefine::scaleWrite(xercesc::DOMElement* element,const G4String& name,const G4ThreeVector& scl) {

   xercesc::DOMElement* scaleElement = newElement("scale");
   element->appendChild(scaleElement);

   if (!name.empty())
      scaleElement->setAttributeNode(newAttribute("name",name));

   G4double x = (fabs(scl.x()-1.0) < kRelativePrecision) ? 1.0 : scl.x();
   G4double y = (fabs(scl.y()-1.0) < kRelativePrecision) ? 1.0 : scl.y();
   G4double z = (fabs(scl.z()-1.0) < kRelativePrecision) ? 1.0 : scl.z();

   scaleElement->setAttributeNode(newAttribute("x",x));
   scaleElement->setAttributeNode(newAttribute("y",y));
   scaleElement->setAttributeNode(newAttribute("z",z));
}

void G4GDMLWriteDefine::rotationWrite(xercesc::DOMElement* element,const G4String& name,const G4ThreeVector& rot) {

   xercesc::DOMElement* rotationElement = newElement("rotation");
   element->appendChild(rotationElement);

   if (!name.empty())
      rotationElement->setAttributeNode(newAttribute("name",name));

   G4double x = (fabs(rot.x()) < kAngularPrecision) ? 0.0 : rot.x();
   G4double y = (fabs(rot.y()) < kAngularPrecision) ? 0.0 : rot.y();
   G4double z = (fabs(rot.z()) < kAngularPrecision) ? 0.0 : rot.z();

   rotationElement->setAttributeNode(newAttribute("x",x/degree));
   rotationElement->setAttributeNode(newAttribute("y",y/degree));
   rotationElement->setAttributeNode(newAttribute("z",z/degree));
   rotationElement->setAttributeNode(newAttribute("unit","deg"));
}

void G4GDMLWriteDefine::positionWrite(xercesc::DOMElement* element,const G4String& name,const G4ThreeVector& pos) {

   xercesc::DOMElement* positionElement = newElement("position");
   element->appendChild(positionElement);

   if (!name.empty())
      positionElement->setAttributeNode(newAttribute("name",name));

   G4double x = (fabs(pos.x()) < kLinearPrecision) ? 0.0 : pos.x();
   G4double y = (fabs(pos.y()) < kLinearPrecision) ? 0.0 : pos.y();
   G4double z = (fabs(pos.z()) < kLinearPrecision) ? 0.0 : pos.z();

   positionElement->setAttributeNode(newAttribute("x",x/mm));
   positionElement->setAttributeNode(newAttribute("y",y/mm));
   positionElement->setAttributeNode(newAttribute("z",z/mm));
   positionElement->setAttributeNode(newAttribute("unit","mm"));
}

void G4GDMLWriteDefine::firstrotationWrite(xercesc::DOMElement* element,const G4ThreeVector& rot) {

   xercesc::DOMElement* rotationElement = newElement("firstrotation");
   element->appendChild(rotationElement);

   G4double x = (fabs(rot.x()) < kAngularPrecision) ? 0.0 : rot.x();
   G4double y = (fabs(rot.y()) < kAngularPrecision) ? 0.0 : rot.y();
   G4double z = (fabs(rot.z()) < kAngularPrecision) ? 0.0 : rot.z();

   rotationElement->setAttributeNode(newAttribute("x",x/degree));
   rotationElement->setAttributeNode(newAttribute("y",y/degree));
   rotationElement->setAttributeNode(newAttribute("z",z/degree));
   rotationElement->setAttributeNode(newAttribute("unit","deg"));
}

void G4GDMLWriteDefine::firstpositionWrite(xercesc::DOMElement* element,const G4ThreeVector& pos) {

   xercesc::DOMElement* positionElement = newElement("firstposition");
   element->appendChild(positionElement);

   G4double x = (fabs(pos.x()) < kLinearPrecision) ? 0.0 : pos.x();
   G4double y = (fabs(pos.y()) < kLinearPrecision) ? 0.0 : pos.y();
   G4double z = (fabs(pos.z()) < kLinearPrecision) ? 0.0 : pos.z();

   positionElement->setAttributeNode(newAttribute("x",x/mm));
   positionElement->setAttributeNode(newAttribute("y",y/mm));
   positionElement->setAttributeNode(newAttribute("z",z/mm));
   positionElement->setAttributeNode(newAttribute("unit","mm"));
}

void G4GDMLWriteDefine::defineWrite(xercesc::DOMElement* element) {

   G4cout << "Writing definitions..." << G4endl;

   defineElement = newElement("define");
   element->appendChild(defineElement);
}

void G4GDMLWriteDefine::AddPosition(const G4String& name,const G4ThreeVector& P) {

   positionWrite(defineElement,name,P);
}
