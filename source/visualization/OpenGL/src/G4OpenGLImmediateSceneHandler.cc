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
// $Id: G4OpenGLImmediateSceneHandler.cc,v 1.17 2005/09/29 14:27:03 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// Andrew Walkden  10th February 1997
// OpenGL immediate scene - draws immediately to buffer
//                           (saving space on server).

#ifdef G4VIS_BUILD_OPENGL_DRIVER

// Included here - problems with HP compiler if not before other includes?
#include "G4NURBS.hh"

// Here follows a special for Mesa, the OpenGL emulator.  Does not affect
// other OpenGL's, as far as I'm aware.   John Allison 18/9/96.
#define CENTERLINE_CLPP  /* CenterLine C++ workaround: */
// Also seems to be required for HP's CC and AIX xlC, at least.

#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLViewerDataStore.hh"
#include "G4OpenGLTransform3D.hh"

#include "G4OpenGLImmediateSceneHandler.hh"

#include "G4LogicalVolume.hh"

G4OpenGLImmediateSceneHandler::G4OpenGLImmediateSceneHandler (G4VGraphicsSystem& system,
						const G4String& name):
G4OpenGLSceneHandler (system, fSceneIdCount++, name)
{}

G4OpenGLImmediateSceneHandler::~G4OpenGLImmediateSceneHandler ()
{}

#include <iomanip>

void G4OpenGLImmediateSceneHandler::BeginPrimitives
(const G4Transform3D& objectTransformation) {
  G4VSceneHandler::BeginPrimitives (objectTransformation);
  glPushMatrix();
  G4OpenGLTransform3D oglt (objectTransformation);

  /*************************** Check matrix.
  const GLdouble* m = oglt.GetGLMatrix ();
  G4cout << "G4OpenGLTransform3D matrix:";
  for (int i = 0; i < 16; i++) {
    if ((i % 4) == 0) G4cout << '\n';
    G4cout << std::setw (15) << m[i];
  }
  G4cout << G4endl;
  *****************************************/

  glMultMatrixd (oglt.GetGLMatrix ());
}

void G4OpenGLImmediateSceneHandler::EndPrimitives () {
  glPopMatrix();
  if (fReadyForTransients) {
    glFlush ();
  }
  G4VSceneHandler::EndPrimitives ();
}

void G4OpenGLImmediateSceneHandler::BeginModeling () {
  G4VSceneHandler::BeginModeling();
}

void G4OpenGLImmediateSceneHandler::EndModeling () {
  G4VSceneHandler::EndModeling ();
}

void G4OpenGLImmediateSceneHandler::ClearTransientStore () {

  G4VSceneHandler::ClearTransientStore ();

  // Make sure screen corresponds to graphical database...
  if (fpViewer) {
    fpViewer -> SetView ();
    fpViewer -> ClearView ();
    fpViewer -> DrawView ();
  }
}

void G4OpenGLImmediateSceneHandler::RequestPrimitives (const G4VSolid& solid) {
  if (fReadyForTransients) {
    // Always draw transient solids, e.g., hits represented as solids.
    // (As we have no control over the order of drawing of transient
    // objects, we cannot do anything about transparent ones, as
    // below, so always draw them.)
    G4VSceneHandler::RequestPrimitives (solid);
  }
  else {
    // For non-transient (run-duration) objects, ensure transparent
    // objects are drawn last.  The problem of
    // blending/transparency/alpha is quite a tricky one - see History
    // of opengl-V07-01-01/2/3.
    // Get vis attributes - pick up defaults if none.
    const G4VisAttributes* pVA =
      fpViewer -> GetApplicableVisAttributes(fpVisAttribs);
    const G4Colour& c = pVA -> GetColour ();
    G4double opacity = c.GetAlpha ();
    if (!fSecondPass) {
      G4bool transparency_enabled =
	G4OpenGLViewerDataStore::GetTransparencyEnabled(fpViewer);
      if (transparency_enabled && opacity < 1.) {
	// On first pass, transparent objects are not drawn, but flag is set...
	fSecondPassRequested = true;
	return;
      }
    }
    // On second pass, opaque objects are not drwan...
    if (fSecondPass && opacity >= 1.) return;
    // Else invoke base class method...
    G4VSceneHandler::RequestPrimitives (solid);
  }
}

G4int G4OpenGLImmediateSceneHandler::fSceneIdCount = 0;

#endif
