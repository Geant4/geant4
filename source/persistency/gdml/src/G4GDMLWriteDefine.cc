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
// $Id: G4GDMLWriteDefine.cc 110108 2018-05-15 11:46:54Z gcosmo $
//
// class G4GDMLWriteDefine Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLWriteDefine.hh"
#include "G4SystemOfUnits.hh"

const G4double G4GDMLWriteDefine::kRelativePrecision = DBL_EPSILON;
const G4double G4GDMLWriteDefine::kAngularPrecision = DBL_EPSILON;
const G4double G4GDMLWriteDefine::kLinearPrecision = DBL_EPSILON;

G4GDMLWriteDefine::G4GDMLWriteDefine()
  : G4GDMLWrite(), defineElement(0)
{
}

G4GDMLWriteDefine::~G4GDMLWriteDefine()
{
}

G4ThreeVector G4GDMLWriteDefine::GetAngles(const G4RotationMatrix& mtx)
{
   G4double x,y,z;
   G4RotationMatrix mat = mtx;
   mat.rectify();   // Rectify matrix from possible roundoff errors

   // Direction of rotation given by left-hand rule; clockwise rotation

   static const G4double kMatrixPrecision = 10E-10;
   const G4double cosb = std::sqrt(mtx.xx()*mtx.xx()+mtx.yx()*mtx.yx());

   if (cosb > kMatrixPrecision)
   {
      x = std::atan2(mtx.zy(),mtx.zz());
      y = std::atan2(-mtx.zx(),cosb);
      z = std::atan2(mtx.yx(),mtx.xx());
   }
   else
   {
      x = std::atan2(-mtx.yz(),mtx.yy());
      y = std::atan2(-mtx.zx(),cosb);
      z = 0.0;
   }

   return G4ThreeVector(x,y,z);
}

void G4GDMLWriteDefine::
Scale_vectorWrite(xercesc::DOMElement* element, const G4String& tag,
                  const G4String& name, const G4ThreeVector& scl)
{
   const G4double x = (std::fabs(scl.x()-1.0) < kRelativePrecision)
                    ? 1.0 : scl.x();
   const G4double y = (std::fabs(scl.y()-1.0) < kRelativePrecision)
                    ? 1.0 : scl.y();
   const G4double z = (std::fabs(scl.z()-1.0) < kRelativePrecision)
                    ? 1.0 : scl.z();

   xercesc::DOMElement* scaleElement = NewElement(tag);
   scaleElement->setAttributeNode(NewAttribute("name",name));
   scaleElement->setAttributeNode(NewAttribute("x",x));
   scaleElement->setAttributeNode(NewAttribute("y",y));
   scaleElement->setAttributeNode(NewAttribute("z",z));
   element->appendChild(scaleElement);
}

void G4GDMLWriteDefine::
Rotation_vectorWrite(xercesc::DOMElement* element, const G4String& tag,
                     const G4String& name, const G4ThreeVector& rot)
{
   const G4double x = (std::fabs(rot.x()) < kAngularPrecision) ? 0.0 : rot.x();
   const G4double y = (std::fabs(rot.y()) < kAngularPrecision) ? 0.0 : rot.y();
   const G4double z = (std::fabs(rot.z()) < kAngularPrecision) ? 0.0 : rot.z();

   xercesc::DOMElement* rotationElement = NewElement(tag);
   rotationElement->setAttributeNode(NewAttribute("name",name));
   rotationElement->setAttributeNode(NewAttribute("x",x/degree));
   rotationElement->setAttributeNode(NewAttribute("y",y/degree));
   rotationElement->setAttributeNode(NewAttribute("z",z/degree));
   rotationElement->setAttributeNode(NewAttribute("unit","deg"));
   element->appendChild(rotationElement);
}

void G4GDMLWriteDefine::
Position_vectorWrite(xercesc::DOMElement* element, const G4String& tag,
                     const G4String& name, const G4ThreeVector& pos)
{
   const G4double x = (std::fabs(pos.x()) < kLinearPrecision) ? 0.0 : pos.x();
   const G4double y = (std::fabs(pos.y()) < kLinearPrecision) ? 0.0 : pos.y();
   const G4double z = (std::fabs(pos.z()) < kLinearPrecision) ? 0.0 : pos.z();

   xercesc::DOMElement* positionElement = NewElement(tag);
   positionElement->setAttributeNode(NewAttribute("name",name));
   positionElement->setAttributeNode(NewAttribute("x",x/mm));
   positionElement->setAttributeNode(NewAttribute("y",y/mm));
   positionElement->setAttributeNode(NewAttribute("z",z/mm));
   positionElement->setAttributeNode(NewAttribute("unit","mm"));
   element->appendChild(positionElement);
}

void G4GDMLWriteDefine::DefineWrite(xercesc::DOMElement* element)
{
#ifdef G4VERBOSE
   G4cout << "G4GDML: Writing definitions..." << G4endl;
#endif
   defineElement = NewElement("define");
   element->appendChild(defineElement);
}
