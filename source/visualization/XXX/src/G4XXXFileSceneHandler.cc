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
//
// 
// John Allison  7th March 2006
// A template for a file-writing graphics driver.
//?? Lines beginning like this require specialisation for your driver.

#include "G4XXXFileSceneHandler.hh"

#include "G4XXXFileViewer.hh"
#include "G4PhysicalVolumeModel.hh"
#include "G4LogicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polyhedron.hh"
#include "G4UnitsTable.hh"

#include <sstream>

G4int G4XXXFileSceneHandler::fSceneIdCount = 0;
// Counter for XXX scene handlers.

G4XXXFileSceneHandler::G4XXXFileSceneHandler(G4VGraphicsSystem& system,
					     const G4String& name):
  G4VSceneHandler(system, fSceneIdCount++, name)
{}

G4XXXFileSceneHandler::~G4XXXFileSceneHandler() {}

#ifdef G4XXXFileDEBUG
// Useful function...
void G4XXXFileSceneHandler::PrintThings() {
  G4cout <<
    "  with transformation "
         << (void*)fpObjectTransformation;
  if (fpModel) {
    G4cout << " from " << fpModel->GetCurrentDescription()
	   << " (tag " << fpModel->GetCurrentTag()
	   << ')';
  } else {
    G4cout << "(not from a model)";
  }
  G4PhysicalVolumeModel* pPVModel =
    dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
  if (pPVModel) {
    G4cout <<
      "\n  current physical volume: "
           << pPVModel->GetCurrentPV()->GetName() <<
      "\n  current logical volume: "
// There might be a problem with the LV pointer if this is a G4LogicalVolumeModel
           << pPVModel->GetCurrentLV()->GetName() <<
      "\n  current depth of geometry tree: "
           << pPVModel->GetCurrentDepth();
  }
  G4cout << G4endl;
}
#endif

// Note: This function overrides G4VSceneHandler::AddSolid(const
// G4Box&).  You may not want to do this, but this is how it's done if
// you do.  Certain other specific solids may be treated this way -
// see G4VSceneHandler.hh.  The simplest possible driver would *not*
// implement these polymorphic functions, with the effect that the
// default versions in G4VSceneHandler are used, which simply call
// G4VSceneHandler::RequestPrimitives to turn the solid into a
// G4Polyhedron usually.
void G4XXXFileSceneHandler::AddSolid(const G4Box& box) {
#ifdef G4XXXFileDEBUG
  G4cout <<
    "G4XXXFileSceneHandler::AddSolid(const G4Box& box) called for "
	 << box.GetName()
	 << G4endl;
#endif
  //?? Process your box...
  std::ostringstream oss;
  oss << "G4Box(" <<
    G4String
    (G4BestUnit
     (G4ThreeVector
      (box.GetXHalfLength(), box.GetYHalfLength(), box.GetZHalfLength()),
      "Length")).strip() << ')';
  dynamic_cast<G4XXXFileViewer*>(fpViewer)->
    GetFileWriter().WriteItem(oss.str());
}

void G4XXXFileSceneHandler::AddPrimitive(const G4Polyline& polyline) {
#ifdef G4XXXFileDEBUG
  G4cout <<
    "G4XXXFileSceneHandler::AddPrimitive(const G4Polyline& polyline) called.\n"
	 << polyline
	 << G4endl;
#endif
  // Get vis attributes - pick up defaults if none.
  //const G4VisAttributes* pVA =
  //  fpViewer -> GetApplicableVisAttributes (polyline.GetVisAttributes ());
  //?? Process polyline.
  std::ostringstream oss;
  oss << polyline;
  dynamic_cast<G4XXXFileViewer*>(fpViewer)->
    GetFileWriter().WriteItem(oss.str());
}

void G4XXXFileSceneHandler::AddPrimitive(const G4Text& text) {
#ifdef G4XXXFileDEBUG
  G4cout <<
    "G4XXXFileSceneHandler::AddPrimitive(const G4Text& text) called.\n"
	 << text
	 << G4endl;
#endif
  // Get text colour - special method since default text colour is
  // determined by the default text vis attributes, which may be
  // specified independent of default vis attributes of other types of
  // visible objects.
  //const G4Colour& c = GetTextColour (text);  // Picks up default if none.
  //?? Process text.
  std::ostringstream oss;
  oss << text;
  dynamic_cast<G4XXXFileViewer*>(fpViewer)->
    GetFileWriter().WriteItem(oss.str());
}

void G4XXXFileSceneHandler::AddPrimitive(const G4Circle& circle) {
#ifdef G4XXXFileDEBUG
  G4cout <<
    "G4XXXFileSceneHandler::AddPrimitive(const G4Circle& circle) called.\n  "
	 << circle
	 << G4endl;
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
#endif
  // Get vis attributes - pick up defaults if none.
  //const G4VisAttributes* pVA =
  //  fpViewer -> GetApplicableVisAttributes (circle.GetVisAttributes ());
  //?? Process circle.
  std::ostringstream oss;
  oss << circle;
  dynamic_cast<G4XXXFileViewer*>(fpViewer)->
    GetFileWriter().WriteItem(oss.str());
}

void G4XXXFileSceneHandler::AddPrimitive(const G4Square& square) {
#ifdef G4XXXFileDEBUG
  G4cout <<
    "G4XXXFileSceneHandler::AddPrimitive(const G4Square& square) called.\n"
	 << square
	 << G4endl;
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
#endif
  // Get vis attributes - pick up defaults if none.
  //const G4VisAttributes* pVA =
  //  fpViewer -> GetApplicableVisAttributes (square.GetVisAttributes ());
  //?? Process square.
  std::ostringstream oss;
  oss << square;
  dynamic_cast<G4XXXFileViewer*>(fpViewer)->
    GetFileWriter().WriteItem(oss.str());
}

void G4XXXFileSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) {
#ifdef G4XXXFileDEBUG
  G4cout <<
"G4XXXFileSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) called.\n"
	 << polyhedron
	 << G4endl;
#endif
  //?? Process polyhedron.
  std::ostringstream oss;
  oss << polyhedron;
  dynamic_cast<G4XXXFileViewer*>(fpViewer)->
    GetFileWriter().WriteItem(oss.str());

  //?? Or... here are some ideas for decomposing into polygons...
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

  // Loop through all the facets...

  // Look at G4OpenGLSceneHandler::AddPrimitive(const G4Polyhedron&)
  // for an example of how to get facets out of a G4Polyhedron,
  // including how to cope with triangles if that's a problem.
}
