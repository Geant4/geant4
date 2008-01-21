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

void G4GDMLWriteSolids::boxWrite(xercesc::DOMElement* solidsElement,const G4Box* const box) {

   xercesc::DOMElement* boxElement = newElement("box");
   solidsElement->appendChild(boxElement);

   boxElement->setAttributeNode(newAttribute("name",box->GetName()));
   boxElement->setAttributeNode(newAttribute("x",2.0*box->GetXHalfLength()));
   boxElement->setAttributeNode(newAttribute("y",2.0*box->GetYHalfLength()));
   boxElement->setAttributeNode(newAttribute("z",2.0*box->GetZHalfLength()));
   boxElement->setAttributeNode(newAttribute("lunit","mm"));
}

void G4GDMLWriteSolids::reflectedSolidWrite(xercesc::DOMElement* solidsElement,const G4ReflectedSolid* const reflectedSolid) {

   xercesc::DOMElement* reflectedSolidElement = newElement("reflectedSolid");
   solidsElement->appendChild(reflectedSolidElement);

   G4Scale3D  scale;
   G4Rotate3D rotation;
   G4Translate3D  translation;

   reflectedSolid->GetTransform3D().getDecomposition(scale,rotation,translation);
   G4ThreeVector angles = getAngles(rotation.getRotation());

   reflectedSolidElement->setAttributeNode(newAttribute("name",reflectedSolid->GetName()));
   reflectedSolidElement->setAttributeNode(newAttribute("solid",reflectedSolid->GetConstituentMovedSolid()->GetName()));
   reflectedSolidElement->setAttributeNode(newAttribute("sx",scale.xx()));
   reflectedSolidElement->setAttributeNode(newAttribute("sy",scale.yy()));
   reflectedSolidElement->setAttributeNode(newAttribute("sz",scale.zz()));
   reflectedSolidElement->setAttributeNode(newAttribute("rx",angles.x()));
   reflectedSolidElement->setAttributeNode(newAttribute("ry",angles.y()));
   reflectedSolidElement->setAttributeNode(newAttribute("rz",angles.z()));
   reflectedSolidElement->setAttributeNode(newAttribute("tx",translation.dx()));
   reflectedSolidElement->setAttributeNode(newAttribute("ty",translation.dy()));
   reflectedSolidElement->setAttributeNode(newAttribute("tz",translation.dz()));
   reflectedSolidElement->setAttributeNode(newAttribute("lunit","mm"));
   reflectedSolidElement->setAttributeNode(newAttribute("aunit","degree"));
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

         addPosition(ref,facet->GetVertex(j));

         NumVertex++;
      }
   }
}

void G4GDMLWriteSolids::solidsWrite(xercesc::DOMElement* element) {

   xercesc::DOMElement* solidsElement = newElement("solids");
   element->appendChild(solidsElement);

   const G4SolidStore* solidList = G4SolidStore::GetInstance();
   const G4int solidCount = solidList->size();

   for (G4int i=0;i<solidCount;i++) {
   
      const G4VSolid* solidPtr = (*solidList)[i];

      if (const G4Box* boxPtr = dynamic_cast<const G4Box*>(solidPtr)) { boxWrite(solidsElement,boxPtr); } else
      if (const G4ReflectedSolid* reflectedSolidPtr = dynamic_cast<const G4ReflectedSolid*>(solidPtr)) { reflectedSolidWrite(solidsElement,reflectedSolidPtr); } else
      if (const G4TessellatedSolid* tessellatedPtr = dynamic_cast<const G4TessellatedSolid*>(solidPtr)) { tessellatedWrite(solidsElement,tessellatedPtr); }
   }
}
