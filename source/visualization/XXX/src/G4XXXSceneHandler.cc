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
// $Id: G4XXXSceneHandler.cc 66870 2013-01-14 23:38:59Z adotti $
//
// 
// John Allison  5th April 2001
// A template for a simplest possible graphics driver.
//?? Lines or sections marked like this require specialisation for your driver.

#include "G4XXXSceneHandler.hh"

#include "G4PhysicalVolumeModel.hh"
#include "G4LogicalVolumeModel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polyhedron.hh"
#include "G4UnitsTable.hh"

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

void G4XXXSceneHandler::AddPrimitive(const G4Polyline&
#ifdef G4XXXDEBUG
 polyline
#endif
) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Polyline& polyline) called.\n"
	 << polyline
	 << G4endl;
  PrintThings();
#endif
  // Get vis attributes - pick up defaults if none.
  //const G4VisAttributes* pVA =
  //  fpViewer -> GetApplicableVisAttributes (polyline.GetVisAttributes ());
  //?? Process polyline.
}

void G4XXXSceneHandler::AddPrimitive(const G4Text&
#ifdef G4XXXDEBUG
 text
#endif
) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Text& text) called.\n"
	 << text
	 << G4endl;
  PrintThings();
#endif
  // Get text colour - special method since default text colour is
  // determined by the default text vis attributes, which may be
  // specified independent of default vis attributes of other types of
  // visible objects.
  //const G4Colour& c = GetTextColour (text);  // Picks up default if none.
  //?? Process text.
}

void G4XXXSceneHandler::AddPrimitive(const G4Circle&
#ifdef G4XXXDEBUG
 circle
#endif
) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Circle& circle) called.\n"
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
  PrintThings();
#endif
  // Get vis attributes - pick up defaults if none.
  //const G4VisAttributes* pVA =
  //  fpViewer -> GetApplicableVisAttributes (circle.GetVisAttributes ());
  //?? Process circle.
}

void G4XXXSceneHandler::AddPrimitive(const G4Square&
#ifdef G4XXXDEBUG
 square
#endif
) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Square& square) called.\n"
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
  PrintThings();
#endif
  // Get vis attributes - pick up defaults if none.
  //const G4VisAttributes* pVA =
  //  fpViewer -> GetApplicableVisAttributes (square.GetVisAttributes ());
  //?? Process square.
}

void G4XXXSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) {
#ifdef G4XXXDEBUG
  G4cout <<
    "G4XXXSceneHandler::AddPrimitive(const G4Polyhedron& polyhedron) called.\n"
	 << polyhedron
	 << G4endl;
  PrintThings();
#endif
  //?? Process polyhedron.  Here are some ideas...
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
