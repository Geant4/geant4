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

G4ThreeVector G4GDMLWriteDefine::getAngles(const G4RotationMatrix& mat) {

   G4double x,y,z;

   G4double cosb = sqrt(mat.xx()*mat.xx()+mat.yx()*mat.yx());

   if (cosb > 16*FLT_EPSILON) {

      x = atan2(mat.zy(),mat.zz());
      y = atan2(-mat.zx(),cosb);
      z = atan2(mat.yx(),mat.xx());
   } else {

      x = atan2(-mat.yz(),mat.yy());
      y = atan2(-mat.zx(),cosb);
      z = 0.0;
   }

   return G4ThreeVector(x,y,z)/CLHEP::degree;
}

void G4GDMLWriteDefine::addPosition(const G4String& name,const G4ThreeVector& P) {

   xercesc::DOMElement* positionElement = newElement("position");
   positionElement->setAttributeNode(newAttribute("name",name));
   positionElement->setAttributeNode(newAttribute("x",P.x()));
   positionElement->setAttributeNode(newAttribute("y",P.y()));
   positionElement->setAttributeNode(newAttribute("z",P.z()));
   positionElement->setAttributeNode(newAttribute("unit","mm"));
   defineElement->appendChild(positionElement);
}

void G4GDMLWriteDefine::positionWrite(xercesc::DOMElement* element,const G4ThreeVector& pos) {

   xercesc::DOMElement* positionElement = newElement("position");
   element->appendChild(positionElement);

   positionElement->setAttributeNode(newAttribute("x",pos.x()));
   positionElement->setAttributeNode(newAttribute("y",pos.y()));
   positionElement->setAttributeNode(newAttribute("z",pos.z()));
   positionElement->setAttributeNode(newAttribute("unit","mm"));
}

void G4GDMLWriteDefine::rotationWrite(xercesc::DOMElement* element,const G4ThreeVector& rot) {

   xercesc::DOMElement* rotationElement = newElement("rotation");
   element->appendChild(rotationElement);

   rotationElement->setAttributeNode(newAttribute("x",rot.x()));
   rotationElement->setAttributeNode(newAttribute("y",rot.y()));
   rotationElement->setAttributeNode(newAttribute("z",rot.z()));
   rotationElement->setAttributeNode(newAttribute("unit","deg"));
}

void G4GDMLWriteDefine::firstpositionWrite(xercesc::DOMElement* element,const G4ThreeVector& pos) {

   xercesc::DOMElement* positionElement = newElement("firstposition");
   element->appendChild(positionElement);

   positionElement->setAttributeNode(newAttribute("x",pos.x()));
   positionElement->setAttributeNode(newAttribute("y",pos.y()));
   positionElement->setAttributeNode(newAttribute("z",pos.z()));
   positionElement->setAttributeNode(newAttribute("unit","mm"));
}

void G4GDMLWriteDefine::firstrotationWrite(xercesc::DOMElement* element,const G4ThreeVector& rot) {

   xercesc::DOMElement* rotationElement = newElement("firstrotation");
   element->appendChild(rotationElement);

   rotationElement->setAttributeNode(newAttribute("x",rot.x()));
   rotationElement->setAttributeNode(newAttribute("y",rot.y()));
   rotationElement->setAttributeNode(newAttribute("z",rot.z()));
   rotationElement->setAttributeNode(newAttribute("unit","deg"));
}

void G4GDMLWriteDefine::scaleWrite(xercesc::DOMElement* element,const G4ThreeVector& scl) {

   xercesc::DOMElement* scaleElement = newElement("scale");
   element->appendChild(scaleElement);

   scaleElement->setAttributeNode(newAttribute("x",scl.x()));
   scaleElement->setAttributeNode(newAttribute("y",scl.y()));
   scaleElement->setAttributeNode(newAttribute("z",scl.z()));
}

void G4GDMLWriteDefine::defineWrite(xercesc::DOMElement* element) {

   G4cout << "Writing definitions..." << G4endl;

   defineElement = newElement("define");
   element->appendChild(defineElement);
}
