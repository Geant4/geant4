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
// $Id: G4GDMLWriteSolids.cc,v 1.48 2008-07-01 08:12:32 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4GDMLWriteSolids Implementation
//
// Original author: Zoltan Torzsok, November 2007
//
// --------------------------------------------------------------------

#include "G4GDMLWriteSolids.hh"

void G4GDMLWriteSolids::booleanWrite(xercesc::DOMElement* solidsElement,const G4BooleanSolid* const boolean) {

   G4String tag("undefined");
   if (dynamic_cast<const G4IntersectionSolid*>(boolean)) { tag = "intersection"; } else
   if (dynamic_cast<const G4SubtractionSolid*>(boolean)) { tag = "subtraction"; } else
   if (dynamic_cast<const G4UnionSolid*>(boolean)) { tag = "union"; }
   
   G4VSolid* firstPtr = const_cast<G4VSolid*>(boolean->GetConstituentSolid(0));
   G4VSolid* secondPtr = const_cast<G4VSolid*>(boolean->GetConstituentSolid(1));
   
   G4ThreeVector firstpos,firstrot,pos,rot;

   if (const G4DisplacedSolid* disp = dynamic_cast<const G4DisplacedSolid*>(firstPtr)) {
   
      firstpos = disp->GetObjectTranslation();
      firstrot = getAngles(disp->GetObjectRotation());
      firstPtr = disp->GetConstituentMovedSolid();
   }

   if (const G4DisplacedSolid* disp = dynamic_cast<const G4DisplacedSolid*>(secondPtr)) {
   
      pos = disp->GetObjectTranslation();
      rot = getAngles(disp->GetObjectRotation());
      secondPtr = disp->GetConstituentMovedSolid();
   }

   AddSolid(firstPtr);   // At first add the constituent solids!
   AddSolid(secondPtr);

   const G4String name = GenerateName(boolean->GetName(),boolean);
   const G4String firstref = GenerateName(firstPtr->GetName(),firstPtr);
   const G4String secondref = GenerateName(secondPtr->GetName(),secondPtr);

   xercesc::DOMElement* booleanElement = newElement(tag);
   booleanElement->setAttributeNode(newAttribute("name",name));
   xercesc::DOMElement* firstElement = newElement("first");
   firstElement->setAttributeNode(newAttribute("ref",firstref));
   booleanElement->appendChild(firstElement);
   xercesc::DOMElement* secondElement = newElement("second");
   secondElement->setAttributeNode(newAttribute("ref",secondref));
   booleanElement->appendChild(secondElement);
   solidsElement->appendChild(booleanElement); // Add the boolean solid AFTER the constituent solids!

   if (fabs(pos.x()) > kLinearPrecision || fabs(pos.y()) > kLinearPrecision || fabs(pos.z()) > kLinearPrecision) positionWrite(booleanElement,name+"position",pos);
   if (fabs(rot.x()) > kAngularPrecision || fabs(rot.y()) > kAngularPrecision || fabs(rot.z()) > kAngularPrecision) rotationWrite(booleanElement,name+"rotation",rot);
   if (fabs(firstpos.x()) > kLinearPrecision || fabs(firstpos.y()) > kLinearPrecision || fabs(firstpos.z()) > kLinearPrecision) firstpositionWrite(booleanElement,name+"firstposition",firstpos);
   if (fabs(firstrot.x()) > kAngularPrecision || fabs(firstrot.y()) > kAngularPrecision || fabs(firstrot.z()) > kAngularPrecision) firstrotationWrite(booleanElement,name+"firstrotation",firstrot);
}

void G4GDMLWriteSolids::boxWrite(xercesc::DOMElement* solidsElement,const G4Box* const box) {

   const G4String name = GenerateName(box->GetName(),box);

   xercesc::DOMElement* boxElement = newElement("box");
   boxElement->setAttributeNode(newAttribute("name",name));
   boxElement->setAttributeNode(newAttribute("x",2.0*box->GetXHalfLength()/mm));
   boxElement->setAttributeNode(newAttribute("y",2.0*box->GetYHalfLength()/mm));
   boxElement->setAttributeNode(newAttribute("z",2.0*box->GetZHalfLength()/mm));
   boxElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(boxElement);
}

void G4GDMLWriteSolids::coneWrite(xercesc::DOMElement* solidsElement,const G4Cons* const cone) {

   const G4String name = GenerateName(cone->GetName(),cone);

   xercesc::DOMElement* coneElement = newElement("cone");
   coneElement->setAttributeNode(newAttribute("name",name));
   coneElement->setAttributeNode(newAttribute("rmin1",cone->GetInnerRadiusMinusZ()/mm));
   coneElement->setAttributeNode(newAttribute("rmax1",cone->GetOuterRadiusMinusZ()/mm));
   coneElement->setAttributeNode(newAttribute("rmin2",cone->GetInnerRadiusPlusZ()/mm));
   coneElement->setAttributeNode(newAttribute("rmax2",cone->GetOuterRadiusPlusZ()/mm));
   coneElement->setAttributeNode(newAttribute("z",2.0*cone->GetZHalfLength()/mm));
   coneElement->setAttributeNode(newAttribute("startphi",cone->GetStartPhiAngle()/degree));
   coneElement->setAttributeNode(newAttribute("deltaphi",cone->GetDeltaPhiAngle()/degree));
   coneElement->setAttributeNode(newAttribute("aunit","deg"));
   coneElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(coneElement);
}

void G4GDMLWriteSolids::ellipsoidWrite(xercesc::DOMElement* solidsElement,const G4Ellipsoid* const ellipsoid) {

   const G4String name = GenerateName(ellipsoid->GetName(),ellipsoid);

   xercesc::DOMElement* ellipsoidElement = newElement("ellipsoid");
   ellipsoidElement->setAttributeNode(newAttribute("name",name));
   ellipsoidElement->setAttributeNode(newAttribute("ax",ellipsoid->GetSemiAxisMax(0)/mm));
   ellipsoidElement->setAttributeNode(newAttribute("by",ellipsoid->GetSemiAxisMax(1)/mm));
   ellipsoidElement->setAttributeNode(newAttribute("cz",ellipsoid->GetSemiAxisMax(2)/mm));
   ellipsoidElement->setAttributeNode(newAttribute("zcut1",ellipsoid->GetZBottomCut()/mm));
   ellipsoidElement->setAttributeNode(newAttribute("zcut2",ellipsoid->GetZTopCut()/mm));
   ellipsoidElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(ellipsoidElement);
}

void G4GDMLWriteSolids::eltubeWrite(xercesc::DOMElement* solidsElement,const G4EllipticalTube* const eltube) {

   const G4String name = GenerateName(eltube->GetName(),eltube);

   xercesc::DOMElement* eltubeElement = newElement("eltube");
   eltubeElement->setAttributeNode(newAttribute("name",name));
   eltubeElement->setAttributeNode(newAttribute("dx",eltube->GetDx()/mm));
   eltubeElement->setAttributeNode(newAttribute("dy",eltube->GetDy()/mm));
   eltubeElement->setAttributeNode(newAttribute("dz",eltube->GetDz()/mm));
   eltubeElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(eltubeElement);
}

void G4GDMLWriteSolids::xtruWrite(xercesc::DOMElement* solidsElement,const G4ExtrudedSolid* const xtru) {

   const G4String name = GenerateName(xtru->GetName(),xtru);

   xercesc::DOMElement* xtruElement = newElement("xtru");
   xtruElement->setAttributeNode(newAttribute("name",name));
   xtruElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(xtruElement);

   const G4int NumVertex = xtru->GetNofVertices();
   
   for (G4int i=0;i<NumVertex;i++) {
   
      xercesc::DOMElement* twoDimVertexElement = newElement("twoDimVertex");
      xtruElement->appendChild(twoDimVertexElement);

      const G4TwoVector vertex = xtru->GetVertex(i);

      twoDimVertexElement->setAttributeNode(newAttribute("x",vertex.x()/mm));
      twoDimVertexElement->setAttributeNode(newAttribute("y",vertex.y()/mm));
   }

   const G4int NumSection = xtru->GetNofZSections();

   for (G4int i=0;i<NumSection;i++) {
   
      xercesc::DOMElement* sectionElement = newElement("section");
      xtruElement->appendChild(sectionElement);

      const G4ExtrudedSolid::ZSection section = xtru->GetZSection(i);

      sectionElement->setAttributeNode(newAttribute("zOrder",i));
      sectionElement->setAttributeNode(newAttribute("zPosition",section.fZ/mm));
      sectionElement->setAttributeNode(newAttribute("xOffset",section.fOffset.x()/mm));
      sectionElement->setAttributeNode(newAttribute("yOffset",section.fOffset.y()/mm));
      sectionElement->setAttributeNode(newAttribute("scalingFactor",section.fScale));
   }
}

void G4GDMLWriteSolids::hypeWrite(xercesc::DOMElement* solidsElement,const G4Hype* const hype) {

   const G4String name = GenerateName(hype->GetName(),hype);

   xercesc::DOMElement* hypeElement = newElement("hype");
   hypeElement->setAttributeNode(newAttribute("name",name));
   hypeElement->setAttributeNode(newAttribute("rmin",hype->GetInnerRadius()/mm));
   hypeElement->setAttributeNode(newAttribute("rmax",hype->GetOuterRadius()/mm));
   hypeElement->setAttributeNode(newAttribute("inst",hype->GetInnerStereo()/degree));
   hypeElement->setAttributeNode(newAttribute("outst",hype->GetOuterStereo()/degree));
   hypeElement->setAttributeNode(newAttribute("z",2.0*hype->GetZHalfLength()/mm));
   hypeElement->setAttributeNode(newAttribute("aunit","deg"));
   hypeElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(hypeElement);
}

void G4GDMLWriteSolids::orbWrite(xercesc::DOMElement* solidsElement,const G4Orb* const orb) {

   const G4String name = GenerateName(orb->GetName(),orb);

   xercesc::DOMElement* orbElement = newElement("orb");
   orbElement->setAttributeNode(newAttribute("name",name));
   orbElement->setAttributeNode(newAttribute("r",orb->GetRadius()/mm));
   orbElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(orbElement);
}

void G4GDMLWriteSolids::paraWrite(xercesc::DOMElement* solidsElement,const G4Para* const para) {

   const G4String name = GenerateName(para->GetName(),para);

   const G4ThreeVector simaxis = para->GetSymAxis();
   const G4double alpha = atan(para->GetTanAlpha());
   const G4double theta = acos(simaxis.z());
   const G4double phi = (simaxis.z() != 1.0) ? (atan(simaxis.y()/simaxis.x())) : (0.0);

   xercesc::DOMElement* paraElement = newElement("para");
   paraElement->setAttributeNode(newAttribute("name",name));
   paraElement->setAttributeNode(newAttribute("x",2.0*para->GetXHalfLength()/mm));
   paraElement->setAttributeNode(newAttribute("y",2.0*para->GetYHalfLength()/mm));
   paraElement->setAttributeNode(newAttribute("z",2.0*para->GetZHalfLength()/mm));
   paraElement->setAttributeNode(newAttribute("alpha",alpha/degree));
   paraElement->setAttributeNode(newAttribute("theta",theta/degree));
   paraElement->setAttributeNode(newAttribute("phi",phi/degree));
   paraElement->setAttributeNode(newAttribute("aunit","deg"));
   paraElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(paraElement);
}

void G4GDMLWriteSolids::polyconeWrite(xercesc::DOMElement* solidsElement,const G4Polycone* const polycone) {

   const G4String name = GenerateName(polycone->GetName(),polycone);

   xercesc::DOMElement* polyconeElement = newElement("polycone");
   polyconeElement->setAttributeNode(newAttribute("name",name));
   polyconeElement->setAttributeNode(newAttribute("startphi",polycone->GetOriginalParameters()->Start_angle/degree));
   polyconeElement->setAttributeNode(newAttribute("deltaphi",polycone->GetOriginalParameters()->Opening_angle/degree));
   polyconeElement->setAttributeNode(newAttribute("aunit","deg"));
   polyconeElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(polyconeElement);

   const size_t num_zplanes = polycone->GetOriginalParameters()->Num_z_planes;
   const G4double* z_array = polycone->GetOriginalParameters()->Z_values;
   const G4double* rmin_array = polycone->GetOriginalParameters()->Rmin;
   const G4double* rmax_array = polycone->GetOriginalParameters()->Rmax;

   for (size_t i=0;i<num_zplanes;i++)
      zplaneWrite(polyconeElement,z_array[i],rmin_array[i],rmax_array[i]);
}

void G4GDMLWriteSolids::polyhedraWrite(xercesc::DOMElement* solidsElement,const G4Polyhedra* const polyhedra) {

   const G4String name = GenerateName(polyhedra->GetName(),polyhedra);

   xercesc::DOMElement* polyhedraElement = newElement("polyhedra");
   polyhedraElement->setAttributeNode(newAttribute("name",name));
   polyhedraElement->setAttributeNode(newAttribute("startphi",polyhedra->GetOriginalParameters()->Start_angle/degree));
   polyhedraElement->setAttributeNode(newAttribute("deltaphi",polyhedra->GetOriginalParameters()->Opening_angle/degree));
   polyhedraElement->setAttributeNode(newAttribute("numsides",polyhedra->GetOriginalParameters()->numSide));
   polyhedraElement->setAttributeNode(newAttribute("aunit","deg"));
   polyhedraElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(polyhedraElement);

   const size_t num_zplanes = polyhedra->GetOriginalParameters()->Num_z_planes;
   const G4double* z_array = polyhedra->GetOriginalParameters()->Z_values;
   const G4double* rmin_array = polyhedra->GetOriginalParameters()->Rmin;
   const G4double* rmax_array = polyhedra->GetOriginalParameters()->Rmax;

   const G4double convertRad = cos(0.5*polyhedra->GetOriginalParameters()->Opening_angle/polyhedra->GetOriginalParameters()->numSide);

   for (size_t i=0;i<num_zplanes;i++)
      zplaneWrite(polyhedraElement,z_array[i],rmin_array[i]*convertRad,rmax_array[i]*convertRad);
}

void G4GDMLWriteSolids::sphereWrite(xercesc::DOMElement* solidsElement,const G4Sphere* const sphere) {

   const G4String name = GenerateName(sphere->GetName(),sphere);

   xercesc::DOMElement* sphereElement = newElement("sphere");
   sphereElement->setAttributeNode(newAttribute("name",name));
   sphereElement->setAttributeNode(newAttribute("rmin",sphere->GetInsideRadius()/mm));
   sphereElement->setAttributeNode(newAttribute("rmax",sphere->GetOuterRadius()/mm));
   sphereElement->setAttributeNode(newAttribute("startphi",sphere->GetStartPhiAngle()/degree));
   sphereElement->setAttributeNode(newAttribute("deltaphi",sphere->GetDeltaPhiAngle()/degree));
   sphereElement->setAttributeNode(newAttribute("starttheta",sphere->GetStartThetaAngle()/degree));
   sphereElement->setAttributeNode(newAttribute("deltatheta",sphere->GetDeltaThetaAngle()/degree));
   sphereElement->setAttributeNode(newAttribute("aunit","deg"));
   sphereElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(sphereElement);
}

void G4GDMLWriteSolids::tessellatedWrite(xercesc::DOMElement* solidsElement,const G4TessellatedSolid* const tessellated) {

   const G4String name = GenerateName(tessellated->GetName(),tessellated);

   xercesc::DOMElement* tessellatedElement = newElement("tessellated");
   tessellatedElement->setAttributeNode(newAttribute("name",name));
   solidsElement->appendChild(tessellatedElement);

   const size_t NumFacets = tessellated->GetNumberOfFacets();
   size_t NumVertex = 0;
   
   for (size_t i=0;i<NumFacets;i++) {
   
      const G4VFacet* facet = tessellated->GetFacet(i);
      const size_t NumVertexPerFacet = facet->GetNumberOfVertices();

      G4String FacetTag;
      
      if (NumVertexPerFacet==3) { FacetTag="triangular"; } else
      if (NumVertexPerFacet==4) { FacetTag="quadrangular"; } else
      G4Exception("G4GDML: ERROR! Facet should contain 3 or 4 vertices!");

      xercesc::DOMElement* facetElement = newElement(FacetTag);
      tessellatedElement->appendChild(facetElement);

      for (size_t j=0;j<NumVertexPerFacet;j++) {
      
         std::stringstream name_stream;
         std::stringstream ref_stream;

         name_stream << "vertex" << (j+1);
	 ref_stream << name << "_vertex" << NumVertex;

         G4String name = name_stream.str();
         G4String ref = ref_stream.str();

         facetElement->setAttributeNode(newAttribute(name,ref));

         AddPosition(ref,facet->GetVertex(j)); // Add position to define section!

         NumVertex++;
      }
   }
}

void G4GDMLWriteSolids::tetWrite(xercesc::DOMElement* solidsElement,const G4Tet* const tet) {

   const G4String name = GenerateName(tet->GetName(),tet);

   std::vector<G4ThreeVector> vertexList = tet->GetVertices();

   xercesc::DOMElement* tetElement = newElement("tet");
   tetElement->setAttributeNode(newAttribute("name",name));
   tetElement->setAttributeNode(newAttribute("vertex1",name+"_vertex1"));
   tetElement->setAttributeNode(newAttribute("vertex2",name+"_vertex2"));
   tetElement->setAttributeNode(newAttribute("vertex3",name+"_vertex3"));
   tetElement->setAttributeNode(newAttribute("vertex4",name+"_vertex4"));
   solidsElement->appendChild(tetElement);

   AddPosition(name+"_vertex1",vertexList[0]);
   AddPosition(name+"_vertex2",vertexList[1]);
   AddPosition(name+"_vertex3",vertexList[2]);
   AddPosition(name+"_vertex4",vertexList[3]);
}

void G4GDMLWriteSolids::torusWrite(xercesc::DOMElement* solidsElement,const G4Torus* const torus) {

   const G4String name = GenerateName(torus->GetName(),torus);

   xercesc::DOMElement* torusElement = newElement("torus");
   torusElement->setAttributeNode(newAttribute("name",name));
   torusElement->setAttributeNode(newAttribute("rmin",torus->GetRmin()/mm));
   torusElement->setAttributeNode(newAttribute("rmax",torus->GetRmax()/mm));
   torusElement->setAttributeNode(newAttribute("rtor",torus->GetRtor()/mm));
   torusElement->setAttributeNode(newAttribute("startphi",torus->GetSPhi()/degree));
   torusElement->setAttributeNode(newAttribute("deltaphi",torus->GetDPhi()/degree));
   torusElement->setAttributeNode(newAttribute("aunit","deg"));
   torusElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(torusElement);
}

void G4GDMLWriteSolids::trapWrite(xercesc::DOMElement* solidsElement,const G4Trap* const trap) {

   const G4String name = GenerateName(trap->GetName(),trap);

   const G4ThreeVector simaxis = trap->GetSymAxis();
   const G4double phi = (simaxis.z() != 1.0) ? (atan(simaxis.y()/simaxis.x())) : (0.0);
   const G4double theta = acos(simaxis.z());
   const G4double alpha1 = atan(trap->GetTanAlpha1());
   const G4double alpha2 = atan(trap->GetTanAlpha2());

   xercesc::DOMElement* trapElement = newElement("trap");
   trapElement->setAttributeNode(newAttribute("name",name));
   trapElement->setAttributeNode(newAttribute("z",2.0*trap->GetZHalfLength()/mm));
   trapElement->setAttributeNode(newAttribute("theta",theta/degree));
   trapElement->setAttributeNode(newAttribute("phi",phi/degree));
   trapElement->setAttributeNode(newAttribute("y1",2.0*trap->GetYHalfLength1()/mm));
   trapElement->setAttributeNode(newAttribute("x1",2.0*trap->GetXHalfLength1()/mm));
   trapElement->setAttributeNode(newAttribute("x2",2.0*trap->GetXHalfLength2()/mm));
   trapElement->setAttributeNode(newAttribute("alpha1",alpha1/degree));
   trapElement->setAttributeNode(newAttribute("y2",2.0*trap->GetYHalfLength2()/mm));
   trapElement->setAttributeNode(newAttribute("x3",2.0*trap->GetXHalfLength3()/mm));
   trapElement->setAttributeNode(newAttribute("x4",2.0*trap->GetXHalfLength4()/mm));
   trapElement->setAttributeNode(newAttribute("alpha2",alpha2/degree));
   trapElement->setAttributeNode(newAttribute("aunit","deg"));
   trapElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(trapElement);
}

void G4GDMLWriteSolids::trdWrite(xercesc::DOMElement* solidsElement,const G4Trd* const trd) {

   const G4String name = GenerateName(trd->GetName(),trd);

   xercesc::DOMElement* trdElement = newElement("trd");
   trdElement->setAttributeNode(newAttribute("name",name));
   trdElement->setAttributeNode(newAttribute("x1",2.0*trd->GetXHalfLength1()/mm));
   trdElement->setAttributeNode(newAttribute("x2",2.0*trd->GetXHalfLength2()/mm));
   trdElement->setAttributeNode(newAttribute("y1",2.0*trd->GetYHalfLength1()/mm));
   trdElement->setAttributeNode(newAttribute("y2",2.0*trd->GetYHalfLength2()/mm));
   trdElement->setAttributeNode(newAttribute("z",2.0*trd->GetZHalfLength()/mm));
   trdElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(trdElement);
}

void G4GDMLWriteSolids::tubeWrite(xercesc::DOMElement* solidsElement,const G4Tubs* const tube) {

   const G4String name = GenerateName(tube->GetName(),tube);

   xercesc::DOMElement* tubeElement = newElement("tube");
   tubeElement->setAttributeNode(newAttribute("name",name));
   tubeElement->setAttributeNode(newAttribute("rmin",tube->GetInnerRadius()/mm));
   tubeElement->setAttributeNode(newAttribute("rmax",tube->GetOuterRadius()/mm));
   tubeElement->setAttributeNode(newAttribute("z",2.0*tube->GetZHalfLength()/mm));
   tubeElement->setAttributeNode(newAttribute("startphi",tube->GetStartPhiAngle()/degree));
   tubeElement->setAttributeNode(newAttribute("deltaphi",tube->GetDeltaPhiAngle()/degree));
   tubeElement->setAttributeNode(newAttribute("aunit","deg"));
   tubeElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(tubeElement);
}

void G4GDMLWriteSolids::twistedboxWrite(xercesc::DOMElement* solidsElement,const G4TwistedBox* const twistedbox) {

   const G4String name = GenerateName(twistedbox->GetName(),twistedbox);

   xercesc::DOMElement* twistedboxElement = newElement("twistedbox");
   twistedboxElement->setAttributeNode(newAttribute("name",name));
   twistedboxElement->setAttributeNode(newAttribute("x",2.0*twistedbox->GetXHalfLength()/mm));
   twistedboxElement->setAttributeNode(newAttribute("y",2.0*twistedbox->GetYHalfLength()/mm));
   twistedboxElement->setAttributeNode(newAttribute("z",2.0*twistedbox->GetZHalfLength()/mm));
   twistedboxElement->setAttributeNode(newAttribute("PhiTwist",twistedbox->GetPhiTwist()/degree));
   twistedboxElement->setAttributeNode(newAttribute("aunit","deg"));
   twistedboxElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(twistedboxElement);
}

void G4GDMLWriteSolids::twistedtrapWrite(xercesc::DOMElement* solidsElement,const G4TwistedTrap* const twistedtrap) {

   const G4String name = GenerateName(twistedtrap->GetName(),twistedtrap);

   xercesc::DOMElement* twistedtrapElement = newElement("twistedtrap");
   twistedtrapElement->setAttributeNode(newAttribute("name",name));
   twistedtrapElement->setAttributeNode(newAttribute("y1",2.0*twistedtrap->GetY1HalfLength()/mm));
   twistedtrapElement->setAttributeNode(newAttribute("x1",2.0*twistedtrap->GetX1HalfLength()/mm));
   twistedtrapElement->setAttributeNode(newAttribute("x2",2.0*twistedtrap->GetX2HalfLength()/mm));
   twistedtrapElement->setAttributeNode(newAttribute("y2",2.0*twistedtrap->GetY2HalfLength()/mm));
   twistedtrapElement->setAttributeNode(newAttribute("x3",2.0*twistedtrap->GetX3HalfLength()/mm));
   twistedtrapElement->setAttributeNode(newAttribute("x4",2.0*twistedtrap->GetX4HalfLength()/mm));
   twistedtrapElement->setAttributeNode(newAttribute("z",2.0*twistedtrap->GetZHalfLength()/mm));
   twistedtrapElement->setAttributeNode(newAttribute("Alph",twistedtrap->GetTiltAngleAlpha()/degree));
   twistedtrapElement->setAttributeNode(newAttribute("theta",twistedtrap->GetPolarAngleTheta()/degree));
   twistedtrapElement->setAttributeNode(newAttribute("phi",twistedtrap->GetAzimuthalAnglePhi()/degree));
   twistedtrapElement->setAttributeNode(newAttribute("PhiTwist",twistedtrap->GetPhiTwist()/degree));
   twistedtrapElement->setAttributeNode(newAttribute("aunit","deg"));
   twistedtrapElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(twistedtrapElement);
}

void G4GDMLWriteSolids::twistedtrdWrite(xercesc::DOMElement* solidsElement,const G4TwistedTrd* const twistedtrd) {

   const G4String name = GenerateName(twistedtrd->GetName(),twistedtrd);

   xercesc::DOMElement* twistedtrdElement = newElement("twistedtrd");
   twistedtrdElement->setAttributeNode(newAttribute("name",name));
   twistedtrdElement->setAttributeNode(newAttribute("x1",2.0*twistedtrd->GetX1HalfLength()/mm));
   twistedtrdElement->setAttributeNode(newAttribute("x2",2.0*twistedtrd->GetX2HalfLength()/mm));
   twistedtrdElement->setAttributeNode(newAttribute("y1",2.0*twistedtrd->GetY1HalfLength()/mm));
   twistedtrdElement->setAttributeNode(newAttribute("y2",2.0*twistedtrd->GetY2HalfLength()/mm));
   twistedtrdElement->setAttributeNode(newAttribute("z",2.0*twistedtrd->GetZHalfLength()/mm));
   twistedtrdElement->setAttributeNode(newAttribute("PhiTwist",twistedtrd->GetPhiTwist()/degree));
   twistedtrdElement->setAttributeNode(newAttribute("aunit","deg"));
   twistedtrdElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(twistedtrdElement);
}

void G4GDMLWriteSolids::twistedtubsWrite(xercesc::DOMElement* solidsElement,const G4TwistedTubs* const twistedtubs) {

   const G4String name = GenerateName(twistedtubs->GetName(),twistedtubs);

   xercesc::DOMElement* twistedtubsElement = newElement("twistedtubs");
   twistedtubsElement->setAttributeNode(newAttribute("name",name));
   twistedtubsElement->setAttributeNode(newAttribute("twistedangle",twistedtubs->GetPhiTwist()/degree));
   twistedtubsElement->setAttributeNode(newAttribute("endinnerrad",twistedtubs->GetInnerRadius()/mm));
   twistedtubsElement->setAttributeNode(newAttribute("endouterrad",twistedtubs->GetOuterRadius()/mm));
   twistedtubsElement->setAttributeNode(newAttribute("zlen",2.0*twistedtubs->GetZHalfLength()/mm));
   twistedtubsElement->setAttributeNode(newAttribute("phi",twistedtubs->GetDPhi()/degree));
   twistedtubsElement->setAttributeNode(newAttribute("aunit","deg"));
   twistedtubsElement->setAttributeNode(newAttribute("lunit","mm"));
   solidsElement->appendChild(twistedtubsElement);
}

void G4GDMLWriteSolids::zplaneWrite(xercesc::DOMElement* element,const G4double& z,const G4double& rmin,const G4double& rmax) {

   xercesc::DOMElement* zplaneElement = newElement("zplane");
   zplaneElement->setAttributeNode(newAttribute("z",z/mm));
   zplaneElement->setAttributeNode(newAttribute("rmin",rmin/mm));
   zplaneElement->setAttributeNode(newAttribute("rmax",rmax/mm));
   element->appendChild(zplaneElement);
}

void G4GDMLWriteSolids::solidsWrite(xercesc::DOMElement* gdmlElement) {

   G4cout << "G4GDML: Writing solids..." << G4endl;

   solidsElement = newElement("solids");
   gdmlElement->appendChild(solidsElement);

   solidList.clear();
}

void G4GDMLWriteSolids::AddSolid(const G4VSolid* const solidPtr) {

   for (size_t i=0;i<solidList.size();i++)   // Check if solid is already in the list!
      if (solidList[i] == solidPtr) return;

   solidList.push_back(solidPtr);

   if (const G4BooleanSolid* const booleanPtr = dynamic_cast<const G4BooleanSolid* const>(solidPtr)) { booleanWrite(solidsElement,booleanPtr); } else
   if (const G4Box* const boxPtr = dynamic_cast<const G4Box* const>(solidPtr)) { boxWrite(solidsElement,boxPtr); } else
   if (const G4Cons* const conePtr = dynamic_cast<const G4Cons* const>(solidPtr)) { coneWrite(solidsElement,conePtr); } else
   if (const G4Ellipsoid* const ellipsoidPtr = dynamic_cast<const G4Ellipsoid* const>(solidPtr)) { ellipsoidWrite(solidsElement,ellipsoidPtr); } else
   if (const G4EllipticalTube* const eltubePtr = dynamic_cast<const G4EllipticalTube* const>(solidPtr)) { eltubeWrite(solidsElement,eltubePtr); } else
   if (const G4ExtrudedSolid* const xtruPtr = dynamic_cast<const G4ExtrudedSolid* const>(solidPtr)) { xtruWrite(solidsElement,xtruPtr); } else
   if (const G4Hype* const hypePtr = dynamic_cast<const G4Hype* const>(solidPtr)) { hypeWrite(solidsElement,hypePtr); } else
   if (const G4Orb* const orbPtr = dynamic_cast<const G4Orb* const>(solidPtr)) { orbWrite(solidsElement,orbPtr); } else
   if (const G4Para* const paraPtr = dynamic_cast<const G4Para* const>(solidPtr)) { paraWrite(solidsElement,paraPtr); } else
   if (const G4Polycone* const polyconePtr = dynamic_cast<const G4Polycone* const>(solidPtr)) { polyconeWrite(solidsElement,polyconePtr); } else
   if (const G4Polyhedra* const polyhedraPtr = dynamic_cast<const G4Polyhedra* const>(solidPtr)) { polyhedraWrite(solidsElement,polyhedraPtr); } else
   if (const G4Sphere* const spherePtr = dynamic_cast<const G4Sphere* const>(solidPtr)) { sphereWrite(solidsElement,spherePtr); } else
   if (const G4TessellatedSolid* const tessellatedPtr = dynamic_cast<const G4TessellatedSolid* const>(solidPtr)) { tessellatedWrite(solidsElement,tessellatedPtr); } else
   if (const G4Tet* const tetPtr = dynamic_cast<const G4Tet* const>(solidPtr)) { tetWrite(solidsElement,tetPtr); } else
   if (const G4Torus* const torusPtr = dynamic_cast<const G4Torus* const>(solidPtr)) { torusWrite(solidsElement,torusPtr); } else
   if (const G4Trap* const trapPtr = dynamic_cast<const G4Trap* const>(solidPtr)) { trapWrite(solidsElement,trapPtr); } else
   if (const G4Trd* const trdPtr = dynamic_cast<const G4Trd* const>(solidPtr)) { trdWrite(solidsElement,trdPtr); } else
   if (const G4Tubs* const tubePtr = dynamic_cast<const G4Tubs* const>(solidPtr)) { tubeWrite(solidsElement,tubePtr); } else
   if (const G4TwistedBox* const twistedboxPtr = dynamic_cast<const G4TwistedBox* const>(solidPtr)) { twistedboxWrite(solidsElement,twistedboxPtr); } else
   if (const G4TwistedTrap* const twistedtrapPtr = dynamic_cast<const G4TwistedTrap* const>(solidPtr)) { twistedtrapWrite(solidsElement,twistedtrapPtr); } else
   if (const G4TwistedTrd* const twistedtrdPtr = dynamic_cast<const G4TwistedTrd* const>(solidPtr)) { twistedtrdWrite(solidsElement,twistedtrdPtr); } else
   if (const G4TwistedTubs* const twistedtubsPtr = dynamic_cast<const G4TwistedTubs* const>(solidPtr)) { twistedtubsWrite(solidsElement,twistedtubsPtr); } else
   G4Exception("G4GDML: ERROR! Unknown solid: "+solidPtr->GetName());
}
