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

#include "G4GDMLWriteParamvol.hh"

void G4GDMLWriteParamvol::box_dimensionsWrite(xercesc::DOMElement* parametersElement,const G4Box* const box) {

   xercesc::DOMElement* box_dimensionsElement = newElement("box_dimensions");
   parametersElement->appendChild(box_dimensionsElement);
   box_dimensionsElement->setAttributeNode(newAttribute("x",2.0*box->GetXHalfLength()));
   box_dimensionsElement->setAttributeNode(newAttribute("y",2.0*box->GetYHalfLength()));
   box_dimensionsElement->setAttributeNode(newAttribute("z",2.0*box->GetZHalfLength()));
   box_dimensionsElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteParamvol::trd_dimensionsWrite(xercesc::DOMElement* parametersElement,const G4Trd* const trd) {

   xercesc::DOMElement* trd_dimensionsElement = newElement("trd_dimensions");
   parametersElement->appendChild(trd_dimensionsElement);
   trd_dimensionsElement->setAttributeNode(newAttribute("name",trd->GetName()));
   trd_dimensionsElement->setAttributeNode(newAttribute("x1",2.0*trd->GetXHalfLength1()));
   trd_dimensionsElement->setAttributeNode(newAttribute("x2",2.0*trd->GetXHalfLength2()));
   trd_dimensionsElement->setAttributeNode(newAttribute("y1",2.0*trd->GetYHalfLength1()));
   trd_dimensionsElement->setAttributeNode(newAttribute("y2",2.0*trd->GetYHalfLength2()));
   trd_dimensionsElement->setAttributeNode(newAttribute("z",2.0*trd->GetZHalfLength()));
   trd_dimensionsElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteParamvol::cone_dimensionsWrite(xercesc::DOMElement* parametersElement,const G4Cons* const cone) {

   xercesc::DOMElement* cone_dimensionsElement = newElement("cone_dimensions");
   parametersElement->appendChild(cone_dimensionsElement);
   cone_dimensionsElement->setAttributeNode(newAttribute("name",cone->GetName()));
   cone_dimensionsElement->setAttributeNode(newAttribute("rmin1",cone->GetInnerRadiusMinusZ()));
   cone_dimensionsElement->setAttributeNode(newAttribute("rmax1",cone->GetOuterRadiusMinusZ()));
   cone_dimensionsElement->setAttributeNode(newAttribute("rmin2",cone->GetInnerRadiusPlusZ()));
   cone_dimensionsElement->setAttributeNode(newAttribute("rmax2",cone->GetOuterRadiusPlusZ()));
   cone_dimensionsElement->setAttributeNode(newAttribute("z",2.0*cone->GetZHalfLength()));
   cone_dimensionsElement->setAttributeNode(newAttribute("startphi",cone->GetStartPhiAngle()/CLHEP::degree));
   cone_dimensionsElement->setAttributeNode(newAttribute("deltaphi",cone->GetDeltaPhiAngle()/CLHEP::degree));
   cone_dimensionsElement->setAttributeNode(newAttribute("aunit","degree"));
   cone_dimensionsElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteParamvol::parametersWrite(xercesc::DOMElement* paramvolElement,const G4VPhysicalVolume* const paramvol,const G4int& index) {

   paramvol->GetParameterisation()->ComputeTransformation(index,const_cast<G4VPhysicalVolume*>(paramvol));
   
   xercesc::DOMElement* parametersElement = newElement("parameters");
   paramvolElement->appendChild(parametersElement);
   rotationWrite(parametersElement,getAngles(paramvol->GetObjectRotationValue()));
   positionWrite(parametersElement,paramvol->GetObjectTranslation());

   G4VSolid* solid = paramvol->GetLogicalVolume()->GetSolid();

   if (G4Box* box = dynamic_cast<G4Box*>(solid)) {
   
      paramvol->GetParameterisation()->ComputeDimensions(*box,index,const_cast<G4VPhysicalVolume*>(paramvol));
      box_dimensionsWrite(parametersElement,box);
   } else
   if (G4Trd* trd = dynamic_cast<G4Trd*>(solid)) {
   
      paramvol->GetParameterisation()->ComputeDimensions(*trd,index,const_cast<G4VPhysicalVolume*>(paramvol));
      trd_dimensionsWrite(parametersElement,trd);
   } else
   if (G4Cons* cone = dynamic_cast<G4Cons*>(solid)) {
   
      paramvol->GetParameterisation()->ComputeDimensions(*cone,index,const_cast<G4VPhysicalVolume*>(paramvol));
      cone_dimensionsWrite(parametersElement,cone);
   }
}

void G4GDMLWriteParamvol::paramvolWrite(xercesc::DOMElement* volumeElement,const G4VPhysicalVolume* const paramvol) {

   xercesc::DOMElement* paramvolElement = newElement("paramvol");
   volumeElement->appendChild(paramvolElement);

   xercesc::DOMElement* volumerefElement = newElement("volumeref");
   paramvolElement->appendChild(volumerefElement);
   volumerefElement->setAttributeNode(newAttribute("ref",paramvol->GetLogicalVolume()->GetName()));

   const G4int parameterCount = paramvol->GetMultiplicity();

   for (G4int i=0;i<parameterCount;i++)
      parametersWrite(paramvolElement,paramvol,i);
}
