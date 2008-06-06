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

void G4GDMLWriteParamvol::trap_dimensionsWrite(xercesc::DOMElement* parametersElement,const G4Trap* const trap) {

   G4ThreeVector simaxis = trap->GetSymAxis();
   G4double phi = (simaxis.z() != 1.0) ? (atan(simaxis.y()/simaxis.x())) : (0.0);
   G4double theta = acos(simaxis.z());
   G4double alpha1 = atan(trap->GetTanAlpha1());
   G4double alpha2 = atan(trap->GetTanAlpha2());

   xercesc::DOMElement* trap_dimensionsElement = newElement("trap");
   parametersElement->appendChild(trap_dimensionsElement);
   trap_dimensionsElement->setAttributeNode(newAttribute("name",trap->GetName()));
   trap_dimensionsElement->setAttributeNode(newAttribute("z",2.0*trap->GetZHalfLength()));
   trap_dimensionsElement->setAttributeNode(newAttribute("theta",theta/CLHEP::degree));
   trap_dimensionsElement->setAttributeNode(newAttribute("phi",phi/CLHEP::degree));
   trap_dimensionsElement->setAttributeNode(newAttribute("y1",2.0*trap->GetYHalfLength1()));
   trap_dimensionsElement->setAttributeNode(newAttribute("x1",2.0*trap->GetXHalfLength1()));
   trap_dimensionsElement->setAttributeNode(newAttribute("x2",2.0*trap->GetXHalfLength2()));
   trap_dimensionsElement->setAttributeNode(newAttribute("alpha1",alpha1/CLHEP::degree));
   trap_dimensionsElement->setAttributeNode(newAttribute("y2",2.0*trap->GetYHalfLength2()));
   trap_dimensionsElement->setAttributeNode(newAttribute("x3",2.0*trap->GetXHalfLength3()));
   trap_dimensionsElement->setAttributeNode(newAttribute("x4",2.0*trap->GetXHalfLength4()));
   trap_dimensionsElement->setAttributeNode(newAttribute("alpha2",alpha2/CLHEP::degree));
   trap_dimensionsElement->setAttributeNode(newAttribute("aunit","degree"));
   trap_dimensionsElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteParamvol::tube_dimensionsWrite(xercesc::DOMElement* parametersElement,const G4Tubs* const tube) {

   xercesc::DOMElement* tube_dimensionsElement = newElement("tube");
   parametersElement->appendChild(tube_dimensionsElement);
   tube_dimensionsElement->setAttributeNode(newAttribute("name",tube->GetName()));
   tube_dimensionsElement->setAttributeNode(newAttribute("rmin",tube->GetInnerRadius()));
   tube_dimensionsElement->setAttributeNode(newAttribute("rmax",tube->GetOuterRadius()));
   tube_dimensionsElement->setAttributeNode(newAttribute("z",2.0*tube->GetZHalfLength()));
   tube_dimensionsElement->setAttributeNode(newAttribute("startphi",tube->GetStartPhiAngle()/CLHEP::degree));
   tube_dimensionsElement->setAttributeNode(newAttribute("deltaphi",tube->GetDeltaPhiAngle()/CLHEP::degree));
   tube_dimensionsElement->setAttributeNode(newAttribute("aunit","degree"));
   tube_dimensionsElement->setAttributeNode(newAttribute("lunit","mm"));
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

void G4GDMLWriteParamvol::sphere_dimensionsWrite(xercesc::DOMElement* parametersElement,const G4Sphere* const sphere) {

   xercesc::DOMElement* sphere_dimensionsElement = newElement("sphere_dimensions");
   parametersElement->appendChild(sphere_dimensionsElement);
   sphere_dimensionsElement->setAttributeNode(newAttribute("name",sphere->GetName()));
   sphere_dimensionsElement->setAttributeNode(newAttribute("rmin",sphere->GetInsideRadius()));
   sphere_dimensionsElement->setAttributeNode(newAttribute("rmax",sphere->GetOuterRadius()));
   sphere_dimensionsElement->setAttributeNode(newAttribute("startphi",sphere->GetStartPhiAngle()/CLHEP::degree));
   sphere_dimensionsElement->setAttributeNode(newAttribute("deltaphi",sphere->GetDeltaPhiAngle()/CLHEP::degree));
   sphere_dimensionsElement->setAttributeNode(newAttribute("starttheta",sphere->GetStartThetaAngle()/CLHEP::degree));
   sphere_dimensionsElement->setAttributeNode(newAttribute("deltatheta",sphere->GetDeltaThetaAngle()/CLHEP::degree));
   sphere_dimensionsElement->setAttributeNode(newAttribute("aunit","degree"));
   sphere_dimensionsElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteParamvol::orb_dimensionsWrite(xercesc::DOMElement* parametersElement,const G4Orb* const orb) {

   xercesc::DOMElement* orb_dimensionsElement = newElement("orb_dimensions");
   parametersElement->appendChild(orb_dimensionsElement);
   orb_dimensionsElement->setAttributeNode(newAttribute("name",orb->GetName()));
   orb_dimensionsElement->setAttributeNode(newAttribute("r",orb->GetRadius()));
   orb_dimensionsElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteParamvol::torus_dimensionsWrite(xercesc::DOMElement* parametersElement,const G4Torus* const torus) {

   xercesc::DOMElement* torus_dimensionsElement = newElement("torus_dimensions");
   parametersElement->appendChild(torus_dimensionsElement);
   torus_dimensionsElement->setAttributeNode(newAttribute("name",torus->GetName()));
   torus_dimensionsElement->setAttributeNode(newAttribute("rmin",torus->GetRmin()));
   torus_dimensionsElement->setAttributeNode(newAttribute("rmax",torus->GetRmax()));
   torus_dimensionsElement->setAttributeNode(newAttribute("rtor",torus->GetRtor()));
   torus_dimensionsElement->setAttributeNode(newAttribute("startphi",torus->GetSPhi()/CLHEP::degree));
   torus_dimensionsElement->setAttributeNode(newAttribute("deltaphi",torus->GetDPhi()/CLHEP::degree));
   torus_dimensionsElement->setAttributeNode(newAttribute("aunit","degree"));
   torus_dimensionsElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteParamvol::para_dimensionsWrite(xercesc::DOMElement* parametersElement,const G4Para* const para) {

   G4ThreeVector simaxis = para->GetSymAxis();
   G4double alpha = atan(para->GetTanAlpha());
   G4double theta = acos(simaxis.z());
   G4double phi = (simaxis.z() != 1.0) ? (atan(simaxis.y()/simaxis.x())) : (0.0);

   xercesc::DOMElement* para_dimensionsElement = newElement("para_dimensions");
   parametersElement->appendChild(para_dimensionsElement);
   para_dimensionsElement->setAttributeNode(newAttribute("name",para->GetName()));
   para_dimensionsElement->setAttributeNode(newAttribute("x",2.0*para->GetXHalfLength()));
   para_dimensionsElement->setAttributeNode(newAttribute("y",2.0*para->GetYHalfLength()));
   para_dimensionsElement->setAttributeNode(newAttribute("z",2.0*para->GetZHalfLength()));
   para_dimensionsElement->setAttributeNode(newAttribute("alpha",alpha/CLHEP::degree));
   para_dimensionsElement->setAttributeNode(newAttribute("theta",theta/CLHEP::degree));
   para_dimensionsElement->setAttributeNode(newAttribute("phi",phi/CLHEP::degree));
   para_dimensionsElement->setAttributeNode(newAttribute("aunit","degree"));
   para_dimensionsElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteParamvol::hype_dimensionsWrite(xercesc::DOMElement* parametersElement,const G4Hype* const hype) {

   xercesc::DOMElement* hype_dimensionsElement = newElement("hype_dimensions");
   parametersElement->appendChild(hype_dimensionsElement);
   hype_dimensionsElement->setAttributeNode(newAttribute("name",hype->GetName()));
   hype_dimensionsElement->setAttributeNode(newAttribute("rmin",hype->GetInnerRadius()));
   hype_dimensionsElement->setAttributeNode(newAttribute("rmax",hype->GetOuterRadius()));
   hype_dimensionsElement->setAttributeNode(newAttribute("inst",hype->GetInnerStereo()/CLHEP::degree));
   hype_dimensionsElement->setAttributeNode(newAttribute("outst",hype->GetOuterStereo()/CLHEP::degree));
   hype_dimensionsElement->setAttributeNode(newAttribute("z",2.0*hype->GetZHalfLength()));
   hype_dimensionsElement->setAttributeNode(newAttribute("aunit","degree"));
   hype_dimensionsElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteParamvol::parametersWrite(xercesc::DOMElement* paramvolElement,const G4VPhysicalVolume* const paramvol,const G4int& index) {

   paramvol->GetParameterisation()->ComputeTransformation(index,const_cast<G4VPhysicalVolume*>(paramvol));
   
   xercesc::DOMElement* parametersElement = newElement("parameters");
   paramvolElement->appendChild(parametersElement);
   rotationWrite(parametersElement,"",getAngles(paramvol->GetObjectRotationValue()));
   positionWrite(parametersElement,"",paramvol->GetObjectTranslation());

   G4VSolid* solid = paramvol->GetLogicalVolume()->GetSolid();

   if (G4Box* box = dynamic_cast<G4Box*>(solid)) {
   
      paramvol->GetParameterisation()->ComputeDimensions(*box,index,const_cast<G4VPhysicalVolume*>(paramvol));
      box_dimensionsWrite(parametersElement,box);
   } else
   if (G4Trd* trd = dynamic_cast<G4Trd*>(solid)) {
   
      paramvol->GetParameterisation()->ComputeDimensions(*trd,index,const_cast<G4VPhysicalVolume*>(paramvol));
      trd_dimensionsWrite(parametersElement,trd);
   } else
   if (G4Trap* trap = dynamic_cast<G4Trap*>(solid)) {
   
      paramvol->GetParameterisation()->ComputeDimensions(*trap,index,const_cast<G4VPhysicalVolume*>(paramvol));
      trap_dimensionsWrite(parametersElement,trap);
   } else
   if (G4Tubs* tube = dynamic_cast<G4Tubs*>(solid)) {
   
      paramvol->GetParameterisation()->ComputeDimensions(*tube,index,const_cast<G4VPhysicalVolume*>(paramvol));
      tube_dimensionsWrite(parametersElement,tube);
   } else
   if (G4Cons* cone = dynamic_cast<G4Cons*>(solid)) {
   
      paramvol->GetParameterisation()->ComputeDimensions(*cone,index,const_cast<G4VPhysicalVolume*>(paramvol));
      cone_dimensionsWrite(parametersElement,cone);
   } else
   if (G4Sphere* sphere = dynamic_cast<G4Sphere*>(solid)) {
   
      paramvol->GetParameterisation()->ComputeDimensions(*sphere,index,const_cast<G4VPhysicalVolume*>(paramvol));
      sphere_dimensionsWrite(parametersElement,sphere);
   } else
   if (G4Orb* orb = dynamic_cast<G4Orb*>(solid)) {
   
      paramvol->GetParameterisation()->ComputeDimensions(*orb,index,const_cast<G4VPhysicalVolume*>(paramvol));
      orb_dimensionsWrite(parametersElement,orb);
   } else
   if (G4Torus* torus = dynamic_cast<G4Torus*>(solid)) {
   
      paramvol->GetParameterisation()->ComputeDimensions(*torus,index,const_cast<G4VPhysicalVolume*>(paramvol));
      torus_dimensionsWrite(parametersElement,torus);
   } else
   if (G4Para* para = dynamic_cast<G4Para*>(solid)) {
   
      paramvol->GetParameterisation()->ComputeDimensions(*para,index,const_cast<G4VPhysicalVolume*>(paramvol));
      para_dimensionsWrite(parametersElement,para);
   } else
   if (G4Hype* hype = dynamic_cast<G4Hype*>(solid)) {
   
      paramvol->GetParameterisation()->ComputeDimensions(*hype,index,const_cast<G4VPhysicalVolume*>(paramvol));
      hype_dimensionsWrite(parametersElement,hype);
   } else
   G4Exception("GDML Writer: ERROR! Solid '"+solid->GetName()+"' can not be used in parameterised volume!");
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
