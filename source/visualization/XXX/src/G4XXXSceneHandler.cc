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
// $Id: G4XXXSceneHandler.cc,v 1.23 2005/06/07 16:46:33 allison Exp $
// GEANT4 tag $Name: geant4-07-01 $
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

G4XXXSceneHandler::G4XXXSceneHandler(G4VGraphicsSystem& system,
					 const G4String& name):
  G4VSceneHandler(system, fSceneIdCount++, name)
{}

G4XXXSceneHandler::~G4XXXSceneHandler() {}

#ifdef G4XXXDEBUG
void G4XXXSceneHandler::PrintThings() {
  G4cout <<
    "  with transformation "
	 << (void*)fpObjectTransformation
	 << " from " << fpModel->GetCurrentDescription()
	 << " (tag " << fpModel->GetCurrentTag()
	 << ')';
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

#include <vector>
std::vector<std::pair<G4VPhysicalVolume*, G4int> > fPVPath;
typedef
std::vector<std::pair<G4VPhysicalVolume*, G4int> >::const_iterator
PVPath_const_iterator;

void G4XXXSceneHandler::PreAddSolid
(const G4Transform3D&, const G4VisAttributes&) {
  using namespace std;
  G4cout <<
    "Current PV/LV/depth: "
	 << fpCurrentPV->GetName() <<
    "/"
	 << fpCurrentLV->GetName() <<
    "/"
	 << fCurrentDepth
	 << G4endl;
  // How to establish a tree (even when some volumes have been
  // culled)...
  static G4int lastDepth = 0;
  static G4LogicalVolume* lastMotherLV = 0;
  G4int copyNo = fpCurrentPV->GetCopyNo();
  if (fCurrentDepth > lastDepth) {
    while (fCurrentDepth > lastDepth++) {
      fPVPath.push_back(make_pair((G4VPhysicalVolume*)0,0));
    }
    fPVPath.back() = make_pair(fpCurrentPV,copyNo);
  } else if (fCurrentDepth == lastDepth) {
    fPVPath.back() = make_pair(fpCurrentPV,copyNo);
  } else {
    while (fCurrentDepth < lastDepth--) {
      fPVPath.pop_back();
    }
    fPVPath.back() = make_pair(fpCurrentPV,copyNo);
  }
  lastDepth = fCurrentDepth;
  lastMotherLV = fpCurrentPV->GetMotherLogical();

  // Debug printing...
  for (PVPath_const_iterator i = fPVPath.begin(); i != fPVPath.end(); ++i) {
    if ((*i).first) {
    G4cout << '/' << (*i).first->GetName();
    } else {
      G4cout << 0;
    }
    G4cout << ':' << (*i).second;
  }
  G4cout << G4endl;

  /***********************************************
  fLVSet.insert(fpCurrentLV);
  G4LogicalVolume* motherLV = fpCurrentPV->GetMotherLogical();
  if (motherLV) {
    if (fLVSet.find(motherLV) == fLVSet.end()) {
      // Mother not previously encountered - must have been culled.
      G4cout << "Mother LV \"" << motherLV->GetName()
	     << "\" not found." << G4endl;
      // Search back up hierarchy...
      do {
	G4LogicalVolume* possibleMotherLV = 0;
	G4LogicalVolumeStore* Store = G4LogicalVolumeStore::GetInstance();
	for ( size_t LV=0; LV < Store->size(); LV++ ) {
	  G4LogicalVolume* aLogVol = (*Store)[LV];
	  if(aLogVol != this) {  // Don't look for it inside itself!
	    for (G4int daughter=0; daughter < aLogVol->GetNoDaughters(); daughter++ ) {
		  if( aLogVol->GetDaughter(daughter)->GetLogicalVolume()==this )
		    { 
		      // aLogVol is the mother !!!
		      //
		      motherLogVol = aLogVol;
		      break;
		    }
		}
	    }
	}
	motherLV = motherLV->GetPhysicalVolume()->GetMotherLogical(); 	
      } while (motherLV && fLVSet.find(motherLV) == fLVSet.end());
      G4cout << "Mother LV \"" << motherLV->GetName()
	     << "\" found." << G4endl;
    }
  }
  ****************************************/

}

void G4XXXSceneHandler::AddSolid(const G4Box& box) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddSolid(const G4Box& box) called for "
	 << box.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddSolid(box);  // Invoke default action.
}

void G4XXXSceneHandler::AddSolid(const G4Cons& cons) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddSolid(const G4Cons& cons) called for "
	 << cons.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddSolid(cons);  // Invoke default action.
}

void G4XXXSceneHandler::AddSolid(const G4Tubs& tubs) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddSolid(const G4Tubs& tubs) called for "
	 << tubs.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddSolid(tubs);  // Invoke default action.
}

void G4XXXSceneHandler::AddSolid(const G4Trd& trd) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddSolid(const G4Trd& trd) called for "
	 << trd.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddSolid(trd);  // Invoke default action.
}

void G4XXXSceneHandler::AddSolid(const G4Trap& trap) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddSolid(const G4Trap& trap) called for "
	 << trap.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddSolid(trap);  // Invoke default action.
}

void G4XXXSceneHandler::AddSolid(const G4Sphere& sphere) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddSolid(const G4Sphere& sphere) called for "
	 << sphere.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddSolid(sphere);  // Invoke default action.
}

void G4XXXSceneHandler::AddSolid(const G4Para& para) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddSolid(const G4Para& para) called for "
	 << para.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddSolid(para);  // Invoke default action.
}

void G4XXXSceneHandler::AddSolid(const G4Torus& torus) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddSolid(const G4Torus& torus) called for "
	 << torus.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddSolid(torus);  // Invoke default action.
}

void G4XXXSceneHandler::AddSolid(const G4Polycone& polycone) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddSolid(const G4Polycone& polycone) called for "
	 << polycone.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddSolid(polycone);  // Invoke default action.
}

void G4XXXSceneHandler::AddSolid(const G4Polyhedra& polyhedra) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddSolid(const G4Polyhedra& polyhedra) called for "
	 << polyhedra.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddSolid(polyhedra);  // Invoke default action.
}

void G4XXXSceneHandler::AddSolid(const G4VSolid& solid) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddSolid(const G4Solid& solid) called for "
	 << solid.GetName()
	 << G4endl;
  PrintThings();
#endif
  G4VSceneHandler::AddSolid(solid);  // Invoke default action.
}

void G4XXXSceneHandler::AddCompound(const G4VTrajectory& traj) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddCompound(const G4VTrajectory& traj) called."
	 << G4endl;
#endif

  G4VSceneHandler::AddCompound(traj);  // Draw trajectory in good old way for now.

  traj.ShowTrajectory();
  G4cout << G4endl;
}

void G4XXXSceneHandler::AddCompound(const G4VHit& hit) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddCompound(const G4VHit& hit) called."
	 << G4endl;
#endif
  G4VSceneHandler::AddCompound(hit);  // Invoke default action.
}

void G4XXXSceneHandler::AddPrimitive(const G4Polyline& polyline) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Polyline& polyline) called:"
    "\n  polyline: " << polyline
	 << G4endl;
  PrintThings();
#endif
  // Get vis attributes - pick up defaults if none.
  //const G4VisAttributes* pVA =
  //  fpViewer -> GetApplicableVisAttributes (polyline.GetVisAttributes ());
}

void G4XXXSceneHandler::AddPrimitive(const G4Text& text) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Text& text) called:"
    "\n  text: " << text.GetText()
	 << G4endl;
  PrintThings();
#endif
  // Get text colour - special method since default text colour is
  // determined by the default text vis attributes, which may be
  // specified independent of default vis attributes of other types of
  // visible objects.
  //const G4Colour& c = GetTextColour (text);  // Picks up default if none.
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
  // Get vis attributes - pick up defaults if none.
  //const G4VisAttributes* pVA =
  //  fpViewer -> GetApplicableVisAttributes (circle.GetVisAttributes ());
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
  // Get vis attributes - pick up defaults if none.
  //const G4VisAttributes* pVA =
  //  fpViewer -> GetApplicableVisAttributes (square.GetVisAttributes ());
}

void G4XXXSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) called."
	 << G4endl;
  PrintThings();
#endif

  //Assume all facets are convex quadrilaterals.
  //Draw each G4Facet individually
  
  //Get colour, etc..
  if (polyhedron.GetNoFacets() == 0) return;

  // Get vis attributes - pick up defaults if none.
  const G4VisAttributes* pVA =
    fpViewer -> GetApplicableVisAttributes (polyhedron.GetVisAttributes ());

  // Get view parameters that the user can force through the vis
  // attributes, thereby over-riding the current view parameter.
  G4ViewParameters::DrawingStyle drawing_style = GetDrawingStyle (pVA);
  G4bool isAuxEdgeVisible = GetAuxEdgeVisible (pVA);
  
  //Get colour, etc..
  //const G4Colour& c = pVA -> GetColour ();
  
  // Initial action depending on drawing style.
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
    G4Point3D vertex;
    G4int edgeFlag;
    G4int edgeCount = 0;
    do {  // loop over edges...
      notLastEdge = polyhedron.GetNextVertex (vertex, edgeFlag);
      // Check to see if edge is visible or not...
      if (isAuxEdgeVisible) {
	edgeFlag = G4int (true);
      }
      if (edgeFlag) {
	// glEdgeFlag (GL_TRUE);
      } else {
	// glEdgeFlag (GL_FALSE);
      }
      // glVertex3d (vertex.x(), 
      //             vertex.y(),
      //             vertex.z());
      edgeCount++;
    } while (notLastEdge);

    // Duplicate last real vertex if necessary to guarantee quadrilateral....
    while (edgeCount < 4) {
      // glEdgeFlag (GL_FALSE);
      // glVertex3d (vertex.x(),
      //             vertex.y(), 
      //             vertex.z());
      edgeCount++;
    }

  } while (notLastFace);  
}

void G4XXXSceneHandler::AddPrimitive(const G4NURBS&) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4NURBS& nurbs) called."
	 << G4endl;
  PrintThings();
#endif
  // Get vis attributes - pick up defaults if none.
  //const G4VisAttributes* pVA =
  //  fpViewer -> GetApplicableVisAttributes (nurbs.GetVisAttributes ());
}

void G4XXXSceneHandler::ClearTransientStore () {
  G4VSceneHandler::ClearTransientStore ();
  // This is typically called after an update and before drawing hits
  // of the next event.  To simulate the clearing of "transients"
  // (hits, etc.) in a system without a graphical database, the
  // detector is redrawn...
  if (fpViewer) {
    fpViewer -> SetView ();
    fpViewer -> ClearView ();
    fpViewer -> DrawView ();
  }
}
