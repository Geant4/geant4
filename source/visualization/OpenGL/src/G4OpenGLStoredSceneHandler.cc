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
// $Id: G4OpenGLStoredSceneHandler.cc,v 1.23 2005/09/29 14:27:03 allison Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// Andrew Walkden  10th February 1997
// OpenGL stored scene - creates OpenGL display lists.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

// Included here - problems with HP compiler if not before other includes?
#include "G4NURBS.hh"

// Here follows a special for Mesa, the OpenGL emulator.  Does not affect
// other OpenGL's, as far as I'm aware.   John Allison 18/9/96.
#define CENTERLINE_CLPP  /* CenterLine C++ workaround: */
// Also seems to be required for HP's CC and AIX xlC, at least.

#include "G4OpenGLStoredSceneHandler.hh"

#include "G4OpenGLViewerDataStore.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"

G4OpenGLStoredSceneHandler::G4OpenGLStoredSceneHandler (G4VGraphicsSystem& system,
					  const G4String& name):
G4OpenGLSceneHandler (system, fSceneIdCount++, name),
fMemoryForDisplayLists (true),
fTopPODL (0)
{}

G4OpenGLStoredSceneHandler::~G4OpenGLStoredSceneHandler ()
{}

void G4OpenGLStoredSceneHandler::BeginPrimitives
(const G4Transform3D& objectTransformation) {
  
  G4VSceneHandler::BeginPrimitives (objectTransformation);

  if (fMemoryForDisplayLists) {
    fDisplayListId = glGenLists (1);
    if (!fDisplayListId) {  // Could pre-allocate?
      G4cout << "********************* WARNING! ********************\n"
	   <<"Unable to allocate any more display lists in OpenGL.\n "
	   << "      Continuing drawing in IMMEDIATE MODE.\n"
	   << "***************************************************" << G4endl;
      fMemoryForDisplayLists = false;
    }
  }
  if (fMemoryForDisplayLists) {
    if (fReadyForTransients) {
      fTODLList.push_back (fDisplayListId);
      fTODLTransformList.push_back (objectTransformation);
      glDrawBuffer (GL_FRONT);
      glPushMatrix();
      G4OpenGLTransform3D oglt (objectTransformation);
      glMultMatrixd (oglt.GetGLMatrix ());
      glNewList (fDisplayListId, GL_COMPILE_AND_EXECUTE);
    }
    else {
      fPODLList.push_back (fDisplayListId);
      fPODLTransformList.push_back (objectTransformation);
      glNewList (fDisplayListId, GL_COMPILE);
    }
  } else {
    glDrawBuffer (GL_FRONT);
    glPushMatrix();
    G4OpenGLTransform3D oglt (objectTransformation);
    glMultMatrixd (oglt.GetGLMatrix ());
  }
}

void G4OpenGLStoredSceneHandler::EndPrimitives () {
  if (fMemoryForDisplayLists) {
    glEndList();
  }
  if (fReadyForTransients || !fMemoryForDisplayLists) {
    glPopMatrix();
    glFlush ();
    glDrawBuffer (GL_BACK);
  }
  G4VSceneHandler::EndPrimitives ();
}

void G4OpenGLStoredSceneHandler::BeginModeling () {
  G4VSceneHandler::BeginModeling();
  ClearStore();  // ...and all that goes with it.
}

void G4OpenGLStoredSceneHandler::EndModeling () {
  // Make a List which calls the other lists.
  fTopPODL = glGenLists (1);
  if (!fTopPODL) {
    G4cout <<
      "ERROR: G4OpenGLStoredSceneHandler::EndModeling: Failure to allocate"
      "  display List for fTopPODL - try OpenGL Immediated mode."
	   << G4endl;
  }
  else {
    glNewList (fTopPODL, GL_COMPILE); {
      for (size_t i = 0; i < fPODLList.size (); i++) {
	glPushMatrix();
	G4OpenGLTransform3D oglt (fPODLTransformList [i]);
	glMultMatrixd (oglt.GetGLMatrix ());
	glCallList (fPODLList[i]);
	glPopMatrix();
      }
    }
    glEndList ();
  }

  G4VSceneHandler::EndModeling ();
}

void G4OpenGLStoredSceneHandler::ClearStore () {

  G4VSceneHandler::ClearStore ();  // Sets need kernel visit, etc.

  size_t i;

  // Delete OpenGL permanent display lists.
  for (i = 0; i < fPODLList.size (); i++) {
    if (fPODLList [i]) {
      glDeleteLists (fPODLList [i], 1);
    } else {
      G4cerr << "Warning : NULL display List in fPODLList." << G4endl;
    }
  }

  if (fTopPODL) glDeleteLists (fTopPODL, 1);
  fTopPODL = 0;

  // Clear other lists, dictionary, etc.
  fPODLList.clear ();
  fPODLTransformList.clear ();
  fSolidMap.clear ();

  // ...and clear transient store...
  for (i = 0; i < fTODLList.size (); i++) {
    if (fTODLList [i]) {
      glDeleteLists (fTODLList [i], 1);
    } else {
      G4cerr << "Warning : NULL display List in fTODLList." << G4endl;
    }
  }
  fTODLList.clear ();
  fTODLTransformList.clear ();
}

void G4OpenGLStoredSceneHandler::ClearTransientStore () {

  G4VSceneHandler::ClearTransientStore ();

  size_t i;

  // Delete OpenGL transient display lists.
  for (i = 0; i < fTODLList.size (); i++) {
    if (fTODLList [i]) {
      glDeleteLists (fTODLList [i], 1);
    } else {
      G4cerr << "Warning : NULL display List in fTODLList." << G4endl;
    }
  }

  // Clear other lists, dictionary, etc.
  fTODLList.clear ();
  fTODLTransformList.clear ();

  // Make sure screen corresponds to graphical database...
  if (fpViewer) {
    fpViewer -> SetView ();
    fpViewer -> ClearView ();
    fpViewer -> DrawView ();
  }
}

void G4OpenGLStoredSceneHandler::RequestPrimitives (const G4VSolid& solid) {
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

    // If a display list already exists for this solid, re-use it if
    // possible.  We could be smarter, and recognise repeated branches
    // of the geometry hierarchy, for example.  But this a;gorithm
    // should be secure, I think...
    const G4VSolid* pSolid = &solid;
    if (fpCurrentPV &&
	// Provided it is not replicated (because if so, the solid's
	// parameters might have been changed)...
	!(fpCurrentPV -> IsReplicated ()) &&
	// ...and if the solid has already been rendered...
	(fSolidMap.find (pSolid) != fSolidMap.end ())) {
      fDisplayListId = fSolidMap [pSolid];
      fPODLList.push_back (fDisplayListId);
      fPODLTransformList.push_back (*fpObjectTransformation);
    }
    else {
      G4VSceneHandler::RequestPrimitives (solid);
      fSolidMap [pSolid] = fDisplayListId;
    }
  }
}

G4int G4OpenGLStoredSceneHandler::fSceneIdCount = 0;

#endif
