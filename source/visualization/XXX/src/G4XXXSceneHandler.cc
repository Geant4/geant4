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
// $Id: G4XXXSceneHandler.cc,v 1.10 2002-11-11 18:26:35 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A base class for a scene handler to dump geometry hierarchy.
// Based on a provisional G4XXXGraphicsScene (was in modeling).

#include "G4XXXSceneHandler.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Torus.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
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
#include "G4VTrajectory.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"

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

#ifdef G4XXXDEBUG
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
#endif

void G4XXXSceneHandler::AddThis(const G4Box& box) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Box& box) called for "
	 << box.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddThis(box);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Cons& cons) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Cons& cons) called for "
	 << cons.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddThis(cons);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Tubs& tubs) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Tubs& tubs) called for "
	 << tubs.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddThis(tubs);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Trd& trd) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Trd& trd) called for "
	 << trd.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddThis(trd);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Trap& trap) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Trap& trap) called for "
	 << trap.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddThis(trap);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Sphere& sphere) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Sphere& sphere) called for "
	 << sphere.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddThis(sphere);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Para& para) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Para& para) called for "
	 << para.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddThis(para);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Torus& torus) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Torus& torus) called for "
	 << torus.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddThis(torus);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Polycone& polycone) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Polycone& polycone) called for "
	 << polycone.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddThis(polycone);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4Polyhedra& polyhedra) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Polyhedra& polyhedra) called for "
	 << polyhedra.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddThis(polyhedra);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4VSolid& solid) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4Solid& solid) called for "
	 << solid.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddThis(solid);  // Invoke default action.
}

void G4XXXSceneHandler::AddThis(const G4VTrajectory& traj) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4VTrajectory& traj) called."
	 << G4endl;
#endif

  G4VSceneHandler::AddThis(traj);  // Draw trajectory in good old way for now.

  traj.ShowTrajectory();
  G4cout << G4endl;
}

void G4XXXSceneHandler::AddThis(const G4VHit& hit) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddThis(const G4VHit& hit) called."
	 << G4endl;
#endif
  G4VSceneHandler::AddThis(hit);  // Invoke default action.
}

void G4XXXSceneHandler::AddPrimitive(const G4Polyline& polyline) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Polyline& polyline) called:"
    "\n  polyline: " << polyline
	 << G4endl;
  PrintThings();
#endif
}

void G4XXXSceneHandler::AddPrimitive(const G4Text& text) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Text& text) called:"
    "\n  text: " << text.GetText()
	 << G4endl;
  PrintThings();
#endif
}

void G4XXXSceneHandler::AddPrimitive(const G4Circle& circle) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Circle& circle) called:\n  ";
  MarkerSizeType sizeType;
  G4double size = GetMarkerSize (circle, sizeType);
  switch (sizeType) {
  default:
  case screen:
    // Draw in screen coordinates.
    G4cout << "screen";
    break;
  case world:
    // Draw in world coordinates.
    G4cout << "world";
    break;
  }
  G4cout << " size: " << size << G4endl;
  PrintThings();
#endif
}

void G4XXXSceneHandler::AddPrimitive(const G4Square& square) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Square& square) called:\n  ";
  MarkerSizeType sizeType;
  G4double size = GetMarkerSize (square, sizeType);
  switch (sizeType) {
  default:
  case screen:
    // Draw in screen coordinates.
    G4cout << "screen";
    break;
  case world:
    // Draw in world coordinates.
    G4cout << "world";
    break;
  }
  G4cout << " size: " << size << G4endl;
  PrintThings();
#endif
}

void G4XXXSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) called."
	 << G4endl;
  PrintThings();
#endif

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
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4NURBS& nurbs) called."
	 << G4endl;
  PrintThings();
#endif
}

void G4XXXSceneHandler::EstablishSpecials
(G4PhysicalVolumeModel& pvModel) {
  pvModel.DefinePointersToWorkingSpace(&fCurrentDepth,
				       &fpCurrentPV,
				       &fpCurrentLV);
}
