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

   if (dynamic_cast<const G4IntersectionSolid* const>(boolean)) { tag = "intersection"; } else
   if (dynamic_cast<const G4SubtractionSolid* const>(boolean)) { tag = "subtraction"; } else
   if (dynamic_cast<const G4UnionSolid* const>(boolean)) { tag = "union"; } else return;

   G4VSolid* firstPtr = const_cast<G4VSolid*>(boolean->GetConstituentSolid(0));
   G4VSolid* secondPtr = const_cast<G4VSolid*>(boolean->GetConstituentSolid(1));

   G4ThreeVector pos;
   G4ThreeVector rot;
   
   if (const G4DisplacedSolid* disp = dynamic_cast<const G4DisplacedSolid*>(secondPtr)) {
   
      pos = disp->GetObjectTranslation();
      rot = getAngles(disp->GetObjectRotation());
   
      secondPtr = disp->GetConstituentMovedSolid();
   }

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
   coneElement->setAttributeNode(newAttribute("startphi",cone->GetStartPhiAngle()));
   coneElement->setAttributeNode(newAttribute("deltaphi",cone->GetDeltaPhiAngle()));
   coneElement->setAttributeNode(newAttribute("lunit","mm"));
   coneElement->setAttributeNode(newAttribute("aunit","degree"));
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

void G4GDMLWriteSolids::tubeWrite(xercesc::DOMElement* solidsElement,const G4Tubs* const tube) {

   xercesc::DOMElement* tubeElement = newElement("tube");
   solidsElement->appendChild(tubeElement);

   tubeElement->setAttributeNode(newAttribute("name",tube->GetName()));
   tubeElement->setAttributeNode(newAttribute("rmin",tube->GetInnerRadius()));
   tubeElement->setAttributeNode(newAttribute("rmax",tube->GetOuterRadius()));
   tubeElement->setAttributeNode(newAttribute("z",2.0*tube->GetZHalfLength()));
   tubeElement->setAttributeNode(newAttribute("startphi",tube->GetStartPhiAngle()));
   tubeElement->setAttributeNode(newAttribute("deltaphi",tube->GetDeltaPhiAngle()));
   tubeElement->setAttributeNode(newAttribute("lunit","mm"));
   tubeElement->setAttributeNode(newAttribute("aunit","degree"));
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

void G4GDMLWriteSolids::solidsWrite(xercesc::DOMElement* element) {

   xercesc::DOMElement* solidsElement = newElement("solids");
   element->appendChild(solidsElement);

   const G4SolidStore* solidList = G4SolidStore::GetInstance();
   const size_t solidCount = solidList->size();

   for (size_t i=0;i<solidCount;i++) {
   
      const G4VSolid* solidPtr = (*solidList)[i];

      if (const G4BooleanSolid* booleanPtr = dynamic_cast<const G4BooleanSolid*>(solidPtr)) { booleanWrite(solidsElement,booleanPtr); } else
      if (const G4Box* boxPtr = dynamic_cast<const G4Box*>(solidPtr)) { boxWrite(solidsElement,boxPtr); } else
      if (const G4Cons* conePtr = dynamic_cast<const G4Cons*>(solidPtr)) { coneWrite(solidsElement,conePtr); } else
      if (const G4ExtrudedSolid* xtruPtr = dynamic_cast<const G4ExtrudedSolid*>(solidPtr)) { xtruWrite(solidsElement,xtruPtr); } else
      if (const G4TessellatedSolid* tessellatedPtr = dynamic_cast<const G4TessellatedSolid*>(solidPtr)) { tessellatedWrite(solidsElement,tessellatedPtr); } else
      if (const G4Tubs* tubePtr = dynamic_cast<const G4Tubs*>(solidPtr)) { tubeWrite(solidsElement,tubePtr); }
   }
}
