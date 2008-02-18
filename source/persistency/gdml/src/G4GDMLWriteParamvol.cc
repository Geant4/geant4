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
