//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4XXXSceneHandler.cc,v 1.5 2001-11-12 18:22:11 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A base class for a scene handler to dump geometry hierarchy.
// Based on a provisional G4XXXGraphicsScene (was in modeling).

#include "G4XXXSceneHandler.hh"

#include "G4VSolid.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4ModelingParameters.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"

G4int G4XXXSceneHandler::fSceneIdCount = 0;
// Counter for XXX scene handlers.

G4int G4XXXSceneHandler::fSceneCount = 0;
// No. of extanct scene handlers.

G4XXXSceneHandler::G4XXXSceneHandler(G4VGraphicsSystem& system,
					 const G4String& name):
  G4VSceneHandler(system, fSceneIdCount++, name),
  fCurrentDepth                 (0),
  fpCurrentPV                   (0),
  fpCurrentLV                   (0)
{
  fSceneCount++;
}

G4XXXSceneHandler::~G4XXXSceneHandler() {}

void G4XXXSceneHandler::PrintThings() {
  G4cout <<
    "  with transformation "
	 << (void*)fpObjectTransformation;
  if (fpCurrentPV) {
    G4cout <<
      "\n  current physical volume: "
	   << fpCurrentPV->GetName() <<
      "\n  current logical volume: "
	   << fpCurrentLV->GetName() <<
      "\n  current depth of geometry tree: "
	   << fCurrentDepth;
  }
  G4cout << G4endl;
}

void G4XXXSceneHandler::AddThis(const G4Box& box) {
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Box& box) called for "
	 << box.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(box);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Cons& cons) {
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Cons& cons) called for "
	 << cons.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(cons);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Tubs& tubs) {
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Tubs& tubs) called for "
	 << tubs.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(tubs);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Trd& trd) {
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Trd& trd) called for "
	 << trd.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(trd);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Trap& trap) {
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Trap& trap) called for "
	 << trap.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(trap);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Sphere& sphere) {
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Sphere& sphere) called for "
	 << sphere.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(sphere);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Para& para) {
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Para& para) called for "
	 << para.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(para);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Torus& torus) {
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Torus& torus) called for "
	 << torus.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(torus);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Polycone& polycone) {
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Polycone& polycone) called for "
	 << polycone.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(polycone);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Polyhedra& polyhedra) {
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Polyhedra& polyhedra) called for "
	 << polyhedra.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(polyhedra);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4VSolid& solid) {
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Solid& solid) called for "
	 << solid.GetName()
	 << G4endl;
  PrintThings();
  G4VSceneHandler::AddThis(solid);  // Invoke default action.
}


void G4XXXSceneHandler::AddPrimitive(const G4Polyline& polyline) {
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Polyline& polyline) called:"
    "\n  polyline: " << polyline
	 << G4endl;
  PrintThings();
}

void G4XXXSceneHandler::AddPrimitive(const G4Text& text) {
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Text& text) called:"
    "\n  text: " << text.GetText()
	 << G4endl;
  PrintThings();
}

void G4XXXSceneHandler::AddPrimitive(const G4Circle& circle) {
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Circle& circle) called:"
    "\n  radius: " << circle.GetWorldRadius()
	 << G4endl;
  PrintThings();
}

void G4XXXSceneHandler::AddPrimitive(const G4Square& square) {
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Square& square) called:"
    "\n  side: " << square.GetWorldRadius()
	 << G4endl;
  PrintThings();
}

void G4XXXSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) {
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) called."
	 << G4endl;
  PrintThings();

  if (polyhedron.GetNoFacets() == 0) return;

  //Assume all facets are convex quadrilaterals.
  //Draw each G4Facet individually
  
  //Get colour, etc..
  // const G4Colour& c = GetColour (polyhedron);
  
  G4ViewParameters::DrawingStyle drawing_style = GetDrawingStyle (polyhedron);
  switch (drawing_style) {
  case (G4ViewParameters::hsr):
    {
      break;
    }
  case (G4ViewParameters::hlr):
    {
      break;
    }
  case (G4ViewParameters::wireframe):
    {
      break;
    }
  default:
    {
      break;
    }     
  }

  //Loop through all the facets...
  G4bool notLastFace;
  G4Normal3D SurfaceUnitNormal;
  do {  // loop over faces...

    //First, find surface normal for the facet...
    notLastFace = polyhedron.GetNextUnitNormal (SurfaceUnitNormal);
        
    //Loop through the four edges of each G4Facet...
    G4bool notLastEdge;
    G4Point3D vertex[4];
    G4int edgeFlag[4], lastEdgeFlag (true);
    edgeFlag[0] = G4int (true);
    G4int edgeCount = 0;
    // glEdgeFlag (GL_TRUE);
    do {  // loop over edges...
      notLastEdge = polyhedron.GetNextVertex (vertex[edgeCount], 
                                              edgeFlag[edgeCount]);
      // Check to see if edge is visible or not...
      if (edgeFlag[edgeCount] != lastEdgeFlag) {
        lastEdgeFlag = edgeFlag[edgeCount];
        if (edgeFlag[edgeCount]) {
          // glEdgeFlag (GL_TRUE);
        } else {
          // glEdgeFlag (GL_FALSE);
        }
      }
      // glVertex3d (vertex[edgeCount].x(), 
      //             vertex[edgeCount].y(),
      //             vertex[edgeCount].z());
      edgeCount++;
    } while (notLastEdge);

    // Duplicate last real vertex if necessary to guarantee quadrilateral....
    while (edgeCount < 4) {
      vertex[edgeCount] = vertex[edgeCount-1];
      // glVertex3d (vertex[edgeCount].x(),
      //             vertex[edgeCount].y(), 
      //             vertex[edgeCount].z());
      edgeCount++;
    }

  } while (notLastFace);  
}

void G4XXXSceneHandler::AddPrimitive(const G4NURBS& nurbs) {
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4NURBS& nurbs) called."
	 << G4endl;
  PrintThings();
}

void G4XXXSceneHandler::EstablishSpecials
(G4PhysicalVolumeModel& pvModel) {
  pvModel.DefinePointersToWorkingSpace(&fCurrentDepth,
				       &fpCurrentPV,
				       &fpCurrentLV);
}
