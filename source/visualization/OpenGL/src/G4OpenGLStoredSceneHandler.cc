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
// $Id: G4OpenGLStoredSceneHandler.cc,v 1.28 2006/06/29 21:19:16 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
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

#include "G4PhysicalVolumeModel.hh"
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

void G4OpenGLStoredSceneHandler::BeginPrimitives2D()
{
  G4VSceneHandler::BeginPrimitives2D();

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
      fTODLTransformList.push_back (G4Transform3D());  // Identity (ignored).
      glDrawBuffer (GL_FRONT);
      glNewList (fDisplayListId, GL_COMPILE_AND_EXECUTE);
    }
    else {
      fPODLList.push_back (fDisplayListId);
      fPODLTransformList.push_back (G4Transform3D());  // Identity (ignored).
      glNewList (fDisplayListId, GL_COMPILE);
    }
  } else {
    glDrawBuffer (GL_FRONT);
  }
  // Push current 3D world matrices and load identity to define screen
  // coordinates...
  glMatrixMode (GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho (-1., 1., -1., 1., -DBL_MAX, DBL_MAX);
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
}

void G4OpenGLStoredSceneHandler::EndPrimitives2D ()
{
  // Pop current 3D world matrices back again...
  glMatrixMode (GL_PROJECTION);
  glPopMatrix();
  glMatrixMode (GL_MODELVIEW);
  glPopMatrix();

  if (fMemoryForDisplayLists) {
    glEndList();
  }
  if (fReadyForTransients || !fMemoryForDisplayLists) {
    glFlush ();
    glDrawBuffer (GL_BACK);
  }
  G4VSceneHandler::EndPrimitives2D ();
}

void G4OpenGLStoredSceneHandler::BeginModeling () {
  G4VSceneHandler::BeginModeling();
  ClearStore();  // ...and all that goes with it.
  /* Debug...
  fDisplayListId = glGenLists (1);
  G4cout << "OGL::fDisplayListId (start): " << fDisplayListId << G4endl;
  */
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
    glNewList (fTopPODL, GL_COMPILE_AND_EXECUTE); {
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

  /* Debug...
  fDisplayListId = glGenLists (1);
  G4cout << "OGL::fDisplayListId (end): " << fDisplayListId << G4endl;
  */
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

void G4OpenGLStoredSceneHandler::RequestPrimitives (const G4VSolid& solid)
{
  if (fReadyForTransients) {
    // Always draw transient solids, e.g., hits represented as solids.
    // (As we have no control over the order of drawing of transient
    // objects, we cannot do anything about transparent ones, as
    // below, so always draw them.)
    G4VSceneHandler::RequestPrimitives (solid);
    return;
  }

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
    G4bool transparency_enabled = true;
    G4OpenGLViewer* pViewer = dynamic_cast<G4OpenGLViewer*>(fpViewer);
    if (pViewer) transparency_enabled = pViewer->transparency_enabled;
    if (transparency_enabled && opacity < 1.) {
      // On first pass, transparent objects are not drawn, but flag is set...
      fSecondPassRequested = true;
      return;
    }
  }

  // On second pass, opaque objects are not drwan...
  if (fSecondPass && opacity >= 1.) return;

  G4PhysicalVolumeModel* pPVModel =
    dynamic_cast<G4PhysicalVolumeModel*>(fpModel);
 
  if (pPVModel) {
    // If part of the geometry hierarchy, i.e., from a
    // G4PhysicalVolumeModel, check if a display list already exists for
    // this solid, re-use it if possible.  We could be smarter, and
    // recognise repeated branches of the geometry hierarchy, for
    // example.  But this algorithm should be secure, I think...
    const G4VSolid* pSolid = &solid;
    EAxis axis = kRho;
    G4VPhysicalVolume* pCurrentPV = pPVModel->GetCurrentPV();
    if (pCurrentPV -> IsReplicated ()) {
      G4int nReplicas;
      G4double width;
      G4double offset;
      G4bool consuming;
      pCurrentPV->GetReplicationData(axis,nReplicas,width,offset,consuming);
    }
    // Provided it is not parametrised (because if so, the
    // solid's parameters might have been changed)...
    if (!(pCurrentPV -> IsParameterised ()) &&
	// Provided it is not replicated radially (because if so, the
	// solid's parameters will have been changed)...
	!(pCurrentPV -> IsReplicated () && axis == kRho) &&
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
    return;
  }

  // Otherwise invoke base class method...
  G4VSceneHandler::RequestPrimitives (solid);
}

G4int G4OpenGLStoredSceneHandler::fSceneIdCount = 0;

#endif
