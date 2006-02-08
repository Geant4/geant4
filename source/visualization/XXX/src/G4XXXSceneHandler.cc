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
// $Id: G4XXXSceneHandler.cc,v 1.27 2006-02-08 15:42:58 allison Exp $
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
#include "G4VisManager.hh"

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

void G4XXXSceneHandler::EstablishSpecials (G4PhysicalVolumeModel& pvModel) {
  pvModel.DefinePointersToWorkingSpace (&fCurrentDepth,
                                        &fpCurrentPV,
                                        &fpCurrentLV,
                                        &fpCurrentMaterial,
					&fDrawnPVPath);
}

void G4XXXSceneHandler::BeginModeling() {
  G4VSceneHandler::BeginModeling();  // Required: see G4VSceneHandler.hh.
}

void G4XXXSceneHandler::EndModeling() {
  fPVNodeStore.clear();
  fDrawnLVStore.clear();
  G4VSceneHandler::EndModeling();  // Required: see G4VSceneHandler.hh.
}

void G4XXXSceneHandler::PreAddSolid
(const G4Transform3D& objectTransformation, const G4VisAttributes& visAttribs)
{
  G4VSceneHandler::PreAddSolid (objectTransformation, visAttribs);

  // fDrawnPVPath is the path of the current drawn (non-culled) volume
  // in terms of drawn (non-culled) ancesters.  Each node is
  // identified by a PVNodeID object, which is a physical volume and
  // copy number.  It is a vector of PVNodeIDs corresponding to the
  // geometry hierarchy actually selected, i.e., not culled.
  typedef G4PhysicalVolumeModel::G4PhysicalVolumeNodeID PVNodeID;
  typedef std::vector<PVNodeID>::iterator PVPath_iterator;
  typedef std::vector<PVNodeID>::reverse_iterator PVPath_reverse_iterator;
  PVPath_reverse_iterator ri;

#ifdef G4XXXDEBUG
  // Current PV:copyNo/LV/depth...
  G4cout << "\nPV:copyNo/LV/depth: "
	 << fpCurrentPV->GetName() << ':' << fpCurrentPV->GetCopyNo() << '/'
	 << fpCurrentLV->GetName() << '/'
	 << fCurrentDepth
	 << G4endl;
  G4cout << "Path of drawn PVs: ";
  for (PVPath_iterator i = fDrawnPVPath.begin();
       i != fDrawnPVPath.end(); ++i) {
    G4cout << '/' << i->GetPhysicalVolume()->GetName()
	   << ':' << i->GetCopyNo();
  }
  G4cout << G4endl;
  // Find mother.  ri points to mother, if any.
  ri = ++fDrawnPVPath.rbegin();
  if (ri != fDrawnPVPath.rend()) {
    // This volume has a mother.
    G4LogicalVolume* drawnMotherLV =
      ri->GetPhysicalVolume()->GetLogicalVolume();
    G4cout << "Mother "
	   << ri->GetPhysicalVolume()->GetName()
	   << ':' << ri->GetCopyNo()
	   << '/' << drawnMotherLV->GetName();
    if (fDrawnLVStore.find(drawnMotherLV) != fDrawnLVStore.end()) {
      // Mother previously encountered.  Add this volume to
      // appropriate node in scene graph tree.
      G4cout << " previously encountered";
    } else {
      // Mother not previously encountered.  Shouldn't happen!
      G4cout << " not previously encountered.  Shouldn't happen!";
    }
  } else {
    // This volume has no mother.  Must be a top level un-culled
    // volume.  Add to root of scene graph tree
    G4cout << "No mother.";
  }
  G4cout << G4endl;
#endif

  // Store ID of current physical volume (not actually used)...
  fPVNodeStore.insert(fDrawnPVPath.back());

  // Actually, it is enough to store the logical volume of current
  // physical volume...
  fDrawnLVStore.insert
    (fDrawnPVPath.back().GetPhysicalVolume()->GetLogicalVolume());

  // Find mother.  ri points to drawn mother, if any.
  ri = ++fDrawnPVPath.rbegin();
  if (ri != fDrawnPVPath.rend()) {
    // This volume has a mother.
    G4LogicalVolume* drawnMotherLV =
      ri->GetPhysicalVolume()->GetLogicalVolume();
    if (fDrawnLVStore.find(drawnMotherLV) != fDrawnLVStore.end()) {
      // Mother previously encountered.  Add this volume to
      // appropriate node in scene graph tree.
      // ...
    } else {
      // Mother not previously encountered.  Shouldn't happen, since
      // G4PhysicalVolumeModel sends volumes as it encounters them,
      // i.e., mothers before daughters, in its descent of the
      // geometry tree.  Error!
      if (G4VisManager::GetInstance()->GetVerbosity() >= G4VisManager::errors){
	G4cout << "ERROR: G4XXXSceneHandler::PreAddSolid: Mother "
	       << ri->GetPhysicalVolume()->GetName()
	       << ':' << ri->GetCopyNo()
	       << " not previously encountered."
	  "\nShouldn't happen!  Please report to visualization coordinator."
	       << G4endl;
      }
      // Continue anyway.  Add to root of scene graph tree.
      // ...
    }
  } else {
    // This volume has no mother.  Must be a top level un-culled
    // volume.  Add to root of scene graph tree.
    // ...
  }
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
  //G4bool isAuxEdgeVisible = GetAuxEdgeVisible (pVA);
  
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

  // Look at G4OpenGLSceneHandler::AddPrimitive(const G4Polyhedron&)
  // for an example of how to get facets out of a G4Polyhedron,
  // including how to cope with triangles if that's a problem.

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
