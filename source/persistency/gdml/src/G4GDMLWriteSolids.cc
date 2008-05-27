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

#include "G4GDMLWriteSolids.hh"

void G4GDMLWriteSolids::booleanWrite(xercesc::DOMElement* solidsElement,const G4BooleanSolid* const boolean) {

   G4String tag;

   if (dynamic_cast<const G4IntersectionSolid*>(boolean)) { tag = "intersection"; } else
   if (dynamic_cast<const G4SubtractionSolid*>(boolean)) { tag = "subtraction"; } else
   if (dynamic_cast<const G4UnionSolid*>(boolean)) { tag = "union"; }
   else { tag = "undefined"; }
   
   G4VSolid* firstPtr = const_cast<G4VSolid*>(boolean->GetConstituentSolid(0));
   G4VSolid* secondPtr = const_cast<G4VSolid*>(boolean->GetConstituentSolid(1));
   
   G4ThreeVector firstpos,firstrot;
   G4ThreeVector pos,rot;
   bool IsFirstDisplaced = false;

   if (const G4DisplacedSolid* disp = dynamic_cast<const G4DisplacedSolid*>(firstPtr)) {
   
      firstpos = disp->GetObjectTranslation();
      firstrot = getAngles(disp->GetObjectRotation());
      firstPtr = disp->GetConstituentMovedSolid();
      IsFirstDisplaced = true;
   }

   if (const G4DisplacedSolid* disp = dynamic_cast<const G4DisplacedSolid*>(secondPtr)) {
   
      pos = disp->GetObjectTranslation();
      rot = getAngles(disp->GetObjectRotation());
      secondPtr = disp->GetConstituentMovedSolid();
   }

   solidsAdd(firstPtr);   // At first add the referenced solids!
   solidsAdd(secondPtr);

   xercesc::DOMElement* booleanElement = newElement(tag);
   solidsElement->appendChild(booleanElement);
   booleanElement->setAttributeNode(newAttribute("name",boolean->GetName()));
   xercesc::DOMElement* firstElement = newElement("first");
   xercesc::DOMElement* secondElement = newElement("second");
   booleanElement->appendChild(firstElement);
   booleanElement->appendChild(secondElement);
   firstElement->setAttributeNode(newAttribute("ref",firstPtr->GetName()));
   secondElement->setAttributeNode(newAttribute("ref",secondPtr->GetName()));

   positionWrite(booleanElement,pos);
   rotationWrite(booleanElement,rot);

   if (IsFirstDisplaced) {
   
      firstpositionWrite(booleanElement,firstpos);
      firstrotationWrite(booleanElement,firstrot);
   }
}

void G4GDMLWriteSolids::boxWrite(xercesc::DOMElement* solidsElement,const G4Box* const box) {

   xercesc::DOMElement* boxElement = newElement("box");
   solidsElement->appendChild(boxElement);
   boxElement->setAttributeNode(newAttribute("name",box->GetName()));
   boxElement->setAttributeNode(newAttribute("x",2.0*box->GetXHalfLength()));
   boxElement->setAttributeNode(newAttribute("y",2.0*box->GetYHalfLength()));
   boxElement->setAttributeNode(newAttribute("z",2.0*box->GetZHalfLength()));
   boxElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteSolids::coneWrite(xercesc::DOMElement* solidsElement,const G4Cons* const cone) {

   xercesc::DOMElement* coneElement = newElement("cone");
   solidsElement->appendChild(coneElement);
   coneElement->setAttributeNode(newAttribute("name",cone->GetName()));
   coneElement->setAttributeNode(newAttribute("rmin1",cone->GetInnerRadiusMinusZ()));
   coneElement->setAttributeNode(newAttribute("rmax1",cone->GetOuterRadiusMinusZ()));
   coneElement->setAttributeNode(newAttribute("rmin2",cone->GetInnerRadiusPlusZ()));
   coneElement->setAttributeNode(newAttribute("rmax2",cone->GetOuterRadiusPlusZ()));
   coneElement->setAttributeNode(newAttribute("z",2.0*cone->GetZHalfLength()));
   coneElement->setAttributeNode(newAttribute("startphi",cone->GetStartPhiAngle()/CLHEP::degree));
   coneElement->setAttributeNode(newAttribute("deltaphi",cone->GetDeltaPhiAngle()/CLHEP::degree));
   coneElement->setAttributeNode(newAttribute("aunit","degree"));
   coneElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteSolids::ellipsoidWrite(xercesc::DOMElement* solidsElement,const G4Ellipsoid* const ellipsoid) {

   xercesc::DOMElement* ellipsoidElement = newElement("ellipsoid");
   solidsElement->appendChild(ellipsoidElement);
   ellipsoidElement->setAttributeNode(newAttribute("name",ellipsoid->GetName()));
   ellipsoidElement->setAttributeNode(newAttribute("ax",ellipsoid->GetSemiAxisMax(0)));
   ellipsoidElement->setAttributeNode(newAttribute("by",ellipsoid->GetSemiAxisMax(1)));
   ellipsoidElement->setAttributeNode(newAttribute("cz",ellipsoid->GetSemiAxisMax(2)));
   ellipsoidElement->setAttributeNode(newAttribute("zcut1",ellipsoid->GetZBottomCut()));
   ellipsoidElement->setAttributeNode(newAttribute("zcut2",ellipsoid->GetZTopCut()));
   ellipsoidElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteSolids::eltubeWrite(xercesc::DOMElement* solidsElement,const G4EllipticalTube* const eltube) {

   xercesc::DOMElement* eltubeElement = newElement("eltube");
   solidsElement->appendChild(eltubeElement);
   eltubeElement->setAttributeNode(newAttribute("name",eltube->GetName()));
   eltubeElement->setAttributeNode(newAttribute("dx",eltube->GetDx()));
   eltubeElement->setAttributeNode(newAttribute("dy",eltube->GetDy()));
   eltubeElement->setAttributeNode(newAttribute("dz",eltube->GetDz()));
   eltubeElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteSolids::xtruWrite(xercesc::DOMElement* solidsElement,const G4ExtrudedSolid* const xtru) {

   xercesc::DOMElement* xtruElement = newElement("xtru");
   solidsElement->appendChild(xtruElement);

   xtruElement->setAttributeNode(newAttribute("name",xtru->GetName()));
   xtruElement->setAttributeNode(newAttribute("lunit","mm"));

   const G4int NumVertex = xtru->GetNofVertices();
   
   for (G4int i=0;i<NumVertex;i++) {
   
      xercesc::DOMElement* twoDimVertexElement = newElement("twoDimVertex");
      xtruElement->appendChild(twoDimVertexElement);

      const G4TwoVector vertex = xtru->GetVertex(i);

      twoDimVertexElement->setAttributeNode(newAttribute("x",vertex.x()));
      twoDimVertexElement->setAttributeNode(newAttribute("y",vertex.y()));
   }

   const G4int NumSection = xtru->GetNofZSections();

   for (G4int i=0;i<NumSection;i++) {
   
      xercesc::DOMElement* sectionElement = newElement("section");
      xtruElement->appendChild(sectionElement);

      const G4ExtrudedSolid::ZSection section = xtru->GetZSection(i);

      sectionElement->setAttributeNode(newAttribute("zOrder",i));
      sectionElement->setAttributeNode(newAttribute("zPosition",section.fZ));
      sectionElement->setAttributeNode(newAttribute("xOffset",section.fOffset.x()));
      sectionElement->setAttributeNode(newAttribute("yOffset",section.fOffset.y()));
      sectionElement->setAttributeNode(newAttribute("scalingFactor",section.fScale));
   }
}

void G4GDMLWriteSolids::hypeWrite(xercesc::DOMElement* solidsElement,const G4Hype* const hype) {

   xercesc::DOMElement* hypeElement = newElement("hype");
   solidsElement->appendChild(hypeElement);
   hypeElement->setAttributeNode(newAttribute("name",hype->GetName()));
   hypeElement->setAttributeNode(newAttribute("rmin",hype->GetInnerRadius()));
   hypeElement->setAttributeNode(newAttribute("rmax",hype->GetOuterRadius()));
   hypeElement->setAttributeNode(newAttribute("inst",hype->GetInnerStereo()/CLHEP::degree));
   hypeElement->setAttributeNode(newAttribute("outst",hype->GetOuterStereo()/CLHEP::degree));
   hypeElement->setAttributeNode(newAttribute("z",2.0*hype->GetZHalfLength()));
   hypeElement->setAttributeNode(newAttribute("aunit","degree"));
   hypeElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteSolids::orbWrite(xercesc::DOMElement* solidsElement,const G4Orb* const orb) {

   xercesc::DOMElement* orbElement = newElement("orb");
   solidsElement->appendChild(orbElement);
   orbElement->setAttributeNode(newAttribute("name",orb->GetName()));
   orbElement->setAttributeNode(newAttribute("r",orb->GetRadius()));
   orbElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteSolids::paraWrite(xercesc::DOMElement* solidsElement,const G4Para* const para) {

   G4ThreeVector simaxis = para->GetSymAxis();
   G4double alpha = atan(para->GetTanAlpha());
   G4double theta = acos(simaxis.z());
   G4double phi = (simaxis.z() != 1.0) ? (atan(simaxis.y()/simaxis.x())) : (0.0);

   xercesc::DOMElement* paraElement = newElement("para");
   solidsElement->appendChild(paraElement);
   paraElement->setAttributeNode(newAttribute("name",para->GetName()));
   paraElement->setAttributeNode(newAttribute("x",2.0*para->GetXHalfLength()));
   paraElement->setAttributeNode(newAttribute("y",2.0*para->GetYHalfLength()));
   paraElement->setAttributeNode(newAttribute("z",2.0*para->GetZHalfLength()));
   paraElement->setAttributeNode(newAttribute("alpha",alpha/CLHEP::degree));
   paraElement->setAttributeNode(newAttribute("theta",theta/CLHEP::degree));
   paraElement->setAttributeNode(newAttribute("phi",phi/CLHEP::degree));
   paraElement->setAttributeNode(newAttribute("aunit","degree"));
   paraElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteSolids::polyconeWrite(xercesc::DOMElement* solidsElement,const G4Polycone* const polycone) {

   xercesc::DOMElement* polyconeElement = newElement("polycone");
   solidsElement->appendChild(polyconeElement);
   polyconeElement->setAttributeNode(newAttribute("name",polycone->GetName()));
   polyconeElement->setAttributeNode(newAttribute("startphi",polycone->GetOriginalParameters()->Start_angle/CLHEP::degree));
   polyconeElement->setAttributeNode(newAttribute("deltaphi",polycone->GetOriginalParameters()->Opening_angle/CLHEP::degree));
   polyconeElement->setAttributeNode(newAttribute("aunit","degree"));
   polyconeElement->setAttributeNode(newAttribute("lunit","mm"));

   const size_t num_zplanes = polycone->GetOriginalParameters()->Num_z_planes;
   const G4double* z_array = polycone->GetOriginalParameters()->Z_values;
   const G4double* rmin_array = polycone->GetOriginalParameters()->Rmin;
   const G4double* rmax_array = polycone->GetOriginalParameters()->Rmax;

   for (size_t i=0;i<num_zplanes;i++)
      zplaneWrite(polyconeElement,z_array[i],rmin_array[i],rmax_array[i]);
}

void G4GDMLWriteSolids::polyhedraWrite(xercesc::DOMElement* solidsElement,const G4Polyhedra* const polyhedra) {

   xercesc::DOMElement* polyhedraElement = newElement("polyhedra");
   solidsElement->appendChild(polyhedraElement);
   polyhedraElement->setAttributeNode(newAttribute("name",polyhedra->GetName()));
   polyhedraElement->setAttributeNode(newAttribute("startphi",polyhedra->GetOriginalParameters()->Start_angle/CLHEP::degree));
   polyhedraElement->setAttributeNode(newAttribute("deltaphi",polyhedra->GetOriginalParameters()->Opening_angle/CLHEP::degree));
   polyhedraElement->setAttributeNode(newAttribute("numsides",polyhedra->GetOriginalParameters()->numSide));
   polyhedraElement->setAttributeNode(newAttribute("aunit","degree"));
   polyhedraElement->setAttributeNode(newAttribute("lunit","mm"));

   const size_t num_zplanes = polyhedra->GetOriginalParameters()->Num_z_planes;
   const G4double* z_array = polyhedra->GetOriginalParameters()->Z_values;
   const G4double* rmin_array = polyhedra->GetOriginalParameters()->Rmin;
   const G4double* rmax_array = polyhedra->GetOriginalParameters()->Rmax;

   G4double convertRad = cos(0.5*polyhedra->GetOriginalParameters()->Opening_angle/polyhedra->GetOriginalParameters()->numSide);

   for (size_t i=0;i<num_zplanes;i++)
      zplaneWrite(polyhedraElement,z_array[i],rmin_array[i]*convertRad,rmax_array[i]*convertRad);
}

void G4GDMLWriteSolids::sphereWrite(xercesc::DOMElement* solidsElement,const G4Sphere* const sphere) {

   xercesc::DOMElement* sphereElement = newElement("sphere");
   solidsElement->appendChild(sphereElement);
   sphereElement->setAttributeNode(newAttribute("name",sphere->GetName()));
   sphereElement->setAttributeNode(newAttribute("rmin",sphere->GetInsideRadius()));
   sphereElement->setAttributeNode(newAttribute("rmax",sphere->GetOuterRadius()));
   sphereElement->setAttributeNode(newAttribute("startphi",sphere->GetStartPhiAngle()/CLHEP::degree));
   sphereElement->setAttributeNode(newAttribute("deltaphi",sphere->GetDeltaPhiAngle()/CLHEP::degree));
   sphereElement->setAttributeNode(newAttribute("starttheta",sphere->GetStartThetaAngle()/CLHEP::degree));
   sphereElement->setAttributeNode(newAttribute("deltatheta",sphere->GetDeltaThetaAngle()/CLHEP::degree));
   sphereElement->setAttributeNode(newAttribute("aunit","degree"));
   sphereElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteSolids::tessellatedWrite(xercesc::DOMElement* solidsElement,const G4TessellatedSolid* const tessellated) {

   xercesc::DOMElement* tessellatedElement = newElement("tessellated");
   solidsElement->appendChild(tessellatedElement);
   tessellatedElement->setAttributeNode(newAttribute("name",tessellated->GetName()));

   const size_t NumFacets = tessellated->GetNumberOfFacets();
   size_t NumVertex = 0;
   
   for (size_t i=0;i<NumFacets;i++) {
   
      const G4VFacet* facet = tessellated->GetFacet(i);
      const size_t NumVertexPerFacet = facet->GetNumberOfVertices();

      G4String FacetTag;
      
      if (NumVertexPerFacet==3) { FacetTag="triangular"; } else
      if (NumVertexPerFacet==4) { FacetTag="quadrangular"; } else
      G4Exception("GDML WRITER: Facet should contain 3 or 4 vertices!");

      xercesc::DOMElement* facetElement = newElement(FacetTag);
      tessellatedElement->appendChild(facetElement);

      for (size_t j=0;j<NumVertexPerFacet;j++) {
      
         std::stringstream name_stream;
         std::stringstream ref_stream;

         name_stream << "vertex" << (j+1);
	 ref_stream << tessellated->GetName() << "_" << "vertex" << NumVertex;

         G4String name = name_stream.str();
         G4String ref = ref_stream.str();

         facetElement->setAttributeNode(newAttribute(name,ref));

         addPosition(ref,facet->GetVertex(j)); // Add position to define section!

         NumVertex++;
      }
   }
}

void G4GDMLWriteSolids::tetWrite(xercesc::DOMElement* solidsElement,const G4Tet* const tet) {

   G4String tetName = tet->GetName();
   std::vector<G4ThreeVector> vertexList = tet->GetVertices();

   xercesc::DOMElement* tetElement = newElement("tet");
   solidsElement->appendChild(tetElement);
   tetElement->setAttributeNode(newAttribute("name",tetName));
   tetElement->setAttributeNode(newAttribute("vertex1",tetName+"_vertex1"));
   tetElement->setAttributeNode(newAttribute("vertex2",tetName+"_vertex2"));
   tetElement->setAttributeNode(newAttribute("vertex3",tetName+"_vertex3"));
   tetElement->setAttributeNode(newAttribute("vertex4",tetName+"_vertex4"));

   addPosition(tetName+"_vertex1",vertexList[0]);
   addPosition(tetName+"_vertex2",vertexList[1]);
   addPosition(tetName+"_vertex3",vertexList[2]);
   addPosition(tetName+"_vertex4",vertexList[3]);
}

void G4GDMLWriteSolids::torusWrite(xercesc::DOMElement* solidsElement,const G4Torus* const torus) {

   xercesc::DOMElement* torusElement = newElement("torus");
   solidsElement->appendChild(torusElement);
   torusElement->setAttributeNode(newAttribute("name",torus->GetName()));
   torusElement->setAttributeNode(newAttribute("rmin",torus->GetRmin()));
   torusElement->setAttributeNode(newAttribute("rmax",torus->GetRmax()));
   torusElement->setAttributeNode(newAttribute("rtor",torus->GetRtor()));
   torusElement->setAttributeNode(newAttribute("startphi",torus->GetSPhi()/CLHEP::degree));
   torusElement->setAttributeNode(newAttribute("deltaphi",torus->GetDPhi()/CLHEP::degree));
   torusElement->setAttributeNode(newAttribute("aunit","degree"));
   torusElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteSolids::trapWrite(xercesc::DOMElement* solidsElement,const G4Trap* const trap) {

   G4ThreeVector simaxis = trap->GetSymAxis();
   G4double phi = (simaxis.z() != 1.0) ? (atan(simaxis.y()/simaxis.x())) : (0.0);
   G4double theta = acos(simaxis.z());
   G4double alpha1 = atan(trap->GetTanAlpha1());
   G4double alpha2 = atan(trap->GetTanAlpha2());

   xercesc::DOMElement* trapElement = newElement("trap");
   solidsElement->appendChild(trapElement);
   trapElement->setAttributeNode(newAttribute("name",trap->GetName()));
   trapElement->setAttributeNode(newAttribute("z",2.0*trap->GetZHalfLength()));
   trapElement->setAttributeNode(newAttribute("theta",theta/CLHEP::degree));
   trapElement->setAttributeNode(newAttribute("phi",phi/CLHEP::degree));
   trapElement->setAttributeNode(newAttribute("y1",2.0*trap->GetYHalfLength1()));
   trapElement->setAttributeNode(newAttribute("x1",2.0*trap->GetXHalfLength1()));
   trapElement->setAttributeNode(newAttribute("x2",2.0*trap->GetXHalfLength2()));
   trapElement->setAttributeNode(newAttribute("alpha1",alpha1/CLHEP::degree));
   trapElement->setAttributeNode(newAttribute("y2",2.0*trap->GetYHalfLength2()));
   trapElement->setAttributeNode(newAttribute("x3",2.0*trap->GetXHalfLength3()));
   trapElement->setAttributeNode(newAttribute("x4",2.0*trap->GetXHalfLength4()));
   trapElement->setAttributeNode(newAttribute("alpha2",alpha2/CLHEP::degree));
   trapElement->setAttributeNode(newAttribute("aunit","degree"));
   trapElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteSolids::trdWrite(xercesc::DOMElement* solidsElement,const G4Trd* const trd) {

   xercesc::DOMElement* trdElement = newElement("trd");
   solidsElement->appendChild(trdElement);
   trdElement->setAttributeNode(newAttribute("name",trd->GetName()));
   trdElement->setAttributeNode(newAttribute("x1",2.0*trd->GetXHalfLength1()));
   trdElement->setAttributeNode(newAttribute("x2",2.0*trd->GetXHalfLength2()));
   trdElement->setAttributeNode(newAttribute("y1",2.0*trd->GetYHalfLength1()));
   trdElement->setAttributeNode(newAttribute("y2",2.0*trd->GetYHalfLength2()));
   trdElement->setAttributeNode(newAttribute("z",2.0*trd->GetZHalfLength()));
   trdElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteSolids::tubeWrite(xercesc::DOMElement* solidsElement,const G4Tubs* const tube) {

   xercesc::DOMElement* tubeElement = newElement("tube");
   solidsElement->appendChild(tubeElement);
   tubeElement->setAttributeNode(newAttribute("name",tube->GetName()));
   tubeElement->setAttributeNode(newAttribute("rmin",tube->GetInnerRadius()));
   tubeElement->setAttributeNode(newAttribute("rmax",tube->GetOuterRadius()));
   tubeElement->setAttributeNode(newAttribute("z",2.0*tube->GetZHalfLength()));
   tubeElement->setAttributeNode(newAttribute("startphi",tube->GetStartPhiAngle()/CLHEP::degree));
   tubeElement->setAttributeNode(newAttribute("deltaphi",tube->GetDeltaPhiAngle()/CLHEP::degree));
   tubeElement->setAttributeNode(newAttribute("aunit","degree"));
   tubeElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteSolids::twistedboxWrite(xercesc::DOMElement* solidsElement,const G4TwistedBox* const twistedbox) {

   xercesc::DOMElement* twistedboxElement = newElement("twistedbox");
   solidsElement->appendChild(twistedboxElement);
   twistedboxElement->setAttributeNode(newAttribute("name",twistedbox->GetName()));
   twistedboxElement->setAttributeNode(newAttribute("x",2.0*twistedbox->GetXHalfLength()));
   twistedboxElement->setAttributeNode(newAttribute("y",2.0*twistedbox->GetYHalfLength()));
   twistedboxElement->setAttributeNode(newAttribute("z",2.0*twistedbox->GetZHalfLength()));
   twistedboxElement->setAttributeNode(newAttribute("PhiTwist",twistedbox->GetPhiTwist()/CLHEP::degree));
   twistedboxElement->setAttributeNode(newAttribute("aunit","degree"));
   twistedboxElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteSolids::twistedtrapWrite(xercesc::DOMElement* solidsElement,const G4TwistedTrap* const twistedtrap) {

   xercesc::DOMElement* twistedtrapElement = newElement("twistedtrap");
   solidsElement->appendChild(twistedtrapElement);
   twistedtrapElement->setAttributeNode(newAttribute("name",twistedtrap->GetName()));
   twistedtrapElement->setAttributeNode(newAttribute("y1",2.0*twistedtrap->GetY1HalfLength()));
   twistedtrapElement->setAttributeNode(newAttribute("x1",2.0*twistedtrap->GetX1HalfLength()));
   twistedtrapElement->setAttributeNode(newAttribute("x2",2.0*twistedtrap->GetX2HalfLength()));
   twistedtrapElement->setAttributeNode(newAttribute("y2",2.0*twistedtrap->GetY2HalfLength()));
   twistedtrapElement->setAttributeNode(newAttribute("x3",2.0*twistedtrap->GetX3HalfLength()));
   twistedtrapElement->setAttributeNode(newAttribute("x4",2.0*twistedtrap->GetX4HalfLength()));
   twistedtrapElement->setAttributeNode(newAttribute("z",2.0*twistedtrap->GetZHalfLength()));
   twistedtrapElement->setAttributeNode(newAttribute("Alph",twistedtrap->GetTiltAngleAlpha()/CLHEP::degree));
   twistedtrapElement->setAttributeNode(newAttribute("theta",twistedtrap->GetPolarAngleTheta()/CLHEP::degree));
   twistedtrapElement->setAttributeNode(newAttribute("phi",twistedtrap->GetAzimuthalAnglePhi()/CLHEP::degree));
   twistedtrapElement->setAttributeNode(newAttribute("PhiTwist",twistedtrap->GetPhiTwist()/CLHEP::degree));
}

void G4GDMLWriteSolids::twistedtrdWrite(xercesc::DOMElement* solidsElement,const G4TwistedTrd* const twistedtrd) {

   xercesc::DOMElement* twistedtrdElement = newElement("twistedtrd");
   solidsElement->appendChild(twistedtrdElement);
   twistedtrdElement->setAttributeNode(newAttribute("name",twistedtrd->GetName()));
   twistedtrdElement->setAttributeNode(newAttribute("x1",2.0*twistedtrd->GetX1HalfLength()));
   twistedtrdElement->setAttributeNode(newAttribute("x2",2.0*twistedtrd->GetX2HalfLength()));
   twistedtrdElement->setAttributeNode(newAttribute("y1",2.0*twistedtrd->GetY1HalfLength()));
   twistedtrdElement->setAttributeNode(newAttribute("y2",2.0*twistedtrd->GetY2HalfLength()));
   twistedtrdElement->setAttributeNode(newAttribute("z",2.0*twistedtrd->GetZHalfLength()));
   twistedtrdElement->setAttributeNode(newAttribute("PhiTwist",twistedtrd->GetPhiTwist()/CLHEP::degree));
   twistedtrdElement->setAttributeNode(newAttribute("aunit","degree"));
   twistedtrdElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteSolids::twistedtubsWrite(xercesc::DOMElement* solidsElement,const G4TwistedTubs* const twistedtubs) {

   xercesc::DOMElement* twistedtubsElement = newElement("twistedtubs");
   solidsElement->appendChild(twistedtubsElement);
   twistedtubsElement->setAttributeNode(newAttribute("name",twistedtubs->GetName()));
   twistedtubsElement->setAttributeNode(newAttribute("twistedangle",twistedtubs->GetPhiTwist()/CLHEP::degree));
   twistedtubsElement->setAttributeNode(newAttribute("endinnerrad",twistedtubs->GetInnerRadius()));
   twistedtubsElement->setAttributeNode(newAttribute("endouterrad",twistedtubs->GetOuterRadius()));
   twistedtubsElement->setAttributeNode(newAttribute("zlen",2.0*twistedtubs->GetZHalfLength()));
   twistedtubsElement->setAttributeNode(newAttribute("phi",twistedtubs->GetDPhi()/CLHEP::degree));
   twistedtubsElement->setAttributeNode(newAttribute("aunit","degree"));
   twistedtubsElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteSolids::zplaneWrite(xercesc::DOMElement* element,const G4double& z,const G4double& rmin,const G4double& rmax) {

   xercesc::DOMElement* zplaneElement = newElement("zplane");
   element->appendChild(zplaneElement);
   zplaneElement->setAttributeNode(newAttribute("z",z));
   zplaneElement->setAttributeNode(newAttribute("rmin",rmin));
   zplaneElement->setAttributeNode(newAttribute("rmax",rmax));
}

void G4GDMLWriteSolids::solidsWrite(xercesc::DOMElement* gdmlElement) {

   G4cout << "Writing solids..." << G4endl;

   solidsElement = newElement("solids");
   gdmlElement->appendChild(solidsElement);
}

void G4GDMLWriteSolids::solidsAdd(const G4VSolid* solidPtr) {

   for (size_t i=0;i<solidList.size();i++)   // Check if solid is already in the list!
      if (solidList[i] == solidPtr) return;

   solidList.push_back(solidPtr);

   if (const G4BooleanSolid* booleanPtr = dynamic_cast<const G4BooleanSolid*>(solidPtr)) { booleanWrite(solidsElement,booleanPtr); } else
   if (const G4Box* boxPtr = dynamic_cast<const G4Box*>(solidPtr)) { boxWrite(solidsElement,boxPtr); } else
   if (const G4Cons* conePtr = dynamic_cast<const G4Cons*>(solidPtr)) { coneWrite(solidsElement,conePtr); } else
   if (const G4Ellipsoid* ellipsoidPtr = dynamic_cast<const G4Ellipsoid*>(solidPtr)) { ellipsoidWrite(solidsElement,ellipsoidPtr); } else
   if (const G4EllipticalTube* eltubePtr = dynamic_cast<const G4EllipticalTube*>(solidPtr)) { eltubeWrite(solidsElement,eltubePtr); } else
   if (const G4ExtrudedSolid* xtruPtr = dynamic_cast<const G4ExtrudedSolid*>(solidPtr)) { xtruWrite(solidsElement,xtruPtr); } else
   if (const G4Hype* hypePtr = dynamic_cast<const G4Hype*>(solidPtr)) { hypeWrite(solidsElement,hypePtr); } else
   if (const G4Orb* orbPtr = dynamic_cast<const G4Orb*>(solidPtr)) { orbWrite(solidsElement,orbPtr); } else
   if (const G4Para* paraPtr = dynamic_cast<const G4Para*>(solidPtr)) { paraWrite(solidsElement,paraPtr); } else
   if (const G4Polycone* polyconePtr = dynamic_cast<const G4Polycone*>(solidPtr)) { polyconeWrite(solidsElement,polyconePtr); } else
   if (const G4Polyhedra* polyhedraPtr = dynamic_cast<const G4Polyhedra*>(solidPtr)) { polyhedraWrite(solidsElement,polyhedraPtr); } else
   if (const G4Sphere* spherePtr = dynamic_cast<const G4Sphere*>(solidPtr)) { sphereWrite(solidsElement,spherePtr); } else
   if (const G4TessellatedSolid* tessellatedPtr = dynamic_cast<const G4TessellatedSolid*>(solidPtr)) { tessellatedWrite(solidsElement,tessellatedPtr); } else
   if (const G4Tet* tetPtr = dynamic_cast<const G4Tet*>(solidPtr)) { tetWrite(solidsElement,tetPtr); } else
   if (const G4Torus* torusPtr = dynamic_cast<const G4Torus*>(solidPtr)) { torusWrite(solidsElement,torusPtr); } else
   if (const G4Trap* trapPtr = dynamic_cast<const G4Trap*>(solidPtr)) { trapWrite(solidsElement,trapPtr); } else
   if (const G4Trd* trdPtr = dynamic_cast<const G4Trd*>(solidPtr)) { trdWrite(solidsElement,trdPtr); } else
   if (const G4Tubs* tubePtr = dynamic_cast<const G4Tubs*>(solidPtr)) { tubeWrite(solidsElement,tubePtr); } else
   if (const G4TwistedBox* twistedboxPtr = dynamic_cast<const G4TwistedBox*>(solidPtr)) { twistedboxWrite(solidsElement,twistedboxPtr); } else
   if (const G4TwistedTrap* twistedtrapPtr = dynamic_cast<const G4TwistedTrap*>(solidPtr)) { twistedtrapWrite(solidsElement,twistedtrapPtr); } else
   if (const G4TwistedTrd* twistedtrdPtr = dynamic_cast<const G4TwistedTrd*>(solidPtr)) { twistedtrdWrite(solidsElement,twistedtrdPtr); } else
   if (const G4TwistedTubs* twistedtubsPtr = dynamic_cast<const G4TwistedTubs*>(solidPtr)) { twistedtubsWrite(solidsElement,twistedtubsPtr); } else
   G4Exception("GDML Writer: ERROR! Unknown solid: "+solidPtr->GetName());
}
