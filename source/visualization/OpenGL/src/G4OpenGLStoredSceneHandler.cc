// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredSceneHandler.cc,v 1.6 1999-12-15 14:54:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLTransform3D.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4Transform3D.hh"
#include "G4Polyline.hh"
#include "G4Text.hh"
#include "G4Circle.hh"
#include "G4Square.hh"
#include "G4Polyhedron.hh"
#include "G4VisAttributes.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ModelingParameters.hh"

#include "G4OpenGLStoredSceneHandler.hh"

G4OpenGLStoredSceneHandler::G4OpenGLStoredSceneHandler (G4VGraphicsSystem& system,
					  const G4String& name):
G4OpenGLSceneHandler (system, fSceneIdCount++, name),
fMemoryForDisplayLists (true),
fTopPODL (0)
{
  fSceneCount++;
}

G4OpenGLStoredSceneHandler::~G4OpenGLStoredSceneHandler ()
{
  fSceneCount--;
}

void G4OpenGLStoredSceneHandler::BeginPrimitives
(const G4Transform3D& objectTransformation) {
  
  G4VSceneHandler::BeginPrimitives (objectTransformation);

  if (fMemoryForDisplayLists) {
    if (!(fDisplayListId = glGenLists (1))) {  // Could pre-allocate?
      G4cout << "********************* WARNING! ********************\n"
	   <<"Unable to allocate any more display lists in OpenGL.\n "
	   << "      Continuing drawing in IMMEDIATE MODE.\n"
	   << "***************************************************" << G4endl;
      fMemoryForDisplayLists = false;
    }
  }
  if (fMemoryForDisplayLists) {
    if (fReadyForTransients) {
      fTODLList.append (fDisplayListId);
      fTODLTransformList.append (objectTransformation);
      glDrawBuffer (GL_FRONT);
      glPushMatrix();
      G4OpenGLTransform3D oglt (objectTransformation);
      glMultMatrixd (oglt.GetGLMatrix ());
      glNewList (fDisplayListId, GL_COMPILE_AND_EXECUTE);
    }
    else {
      fPODLList.append (fDisplayListId);
      fPODLTransformList.append (objectTransformation);
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

void G4OpenGLStoredSceneHandler::ClearStore () {

  G4VSceneHandler::ClearStore ();  // Sets need kernel visit, etc.

  int i;

  // Delete OpenGL display lists.
  for (i = 0; i < fPODLList.entries (); i++) {
    if (fPODLList (i)) {
      glDeleteLists (fPODLList (i), 1);
    } else {
      G4cerr << "Warning : NULL display List in fPODLList." << G4endl;
    }
  }
  for (i = 0; i < fTODLList.entries (); i++) {
    if (fTODLList (i)) {
      glDeleteLists (fTODLList (i), 1);
    } else {
      G4cerr << "Warning : NULL display List in fTODLList." << G4endl;
    }
  }

  if (fTopPODL) glDeleteLists (fTopPODL, 1);
  fTopPODL = 0;

  fMemoryForDisplayLists = glIsList (fTopPODL = glGenLists (1));

  // Clear other lists, dictionary, etc.
  fPODLList.clear ();
  fTODLList.clear ();
  fPODLTransformList.clear ();
  fTODLTransformList.clear ();
  fSolidMap.clear ();
}

void G4OpenGLStoredSceneHandler::BeginModeling () {

  if (fpViewer -> GetViewParameters ().GetDrawingStyle() == G4ViewParameters::hlr) {
    initialize_hlr = true;
  }
  G4VSceneHandler::BeginModeling();
}

void G4OpenGLStoredSceneHandler::EndModeling () {
  // Make a List which calls the other lists.

//  if (!(fTopPODL = glGenLists (1))) {
//    G4Exception ("Unable to allocate display List for fTopPODL in G4OpenGLStoredSceneHandler");
//  }

  glNewList (fTopPODL, GL_COMPILE); {
    for (int i = 0; i < fPODLList.entries (); i++) {
      glPushMatrix();
      G4OpenGLTransform3D oglt (fPODLTransformList (i));
      glMultMatrixd (oglt.GetGLMatrix ());
      glCallList (fPODLList(i));
      glPopMatrix();
    }
  }
  glEndList ();
  G4VSceneHandler::EndModeling ();

  if (fpViewer -> GetViewParameters ().GetDrawingStyle() == G4ViewParameters::hlr) {
    initialize_hlr = false;
    //    glDisable (GL_POLYGON_OFFSET_FILL);
  }
}

void G4OpenGLStoredSceneHandler::RequestPrimitives (const G4VSolid& solid) {
  if (fReadyForTransients ||
      GetModel () -> GetModelingParameters () -> GetRepStyle () ==
      G4ModelingParameters::hierarchy) {
    G4VSceneHandler::RequestPrimitives (solid);
  }
  else {
    // Stop-gap solution for display List re-use.  A proper
    // implementation would use geometry hierarchy.
    const G4VSolid* pSolid = &solid;
    if (!(fpCurrentPV -> IsReplicated ()) &&
	(fSolidMap.find (pSolid) != fSolidMap.end ())) {
      fDisplayListId = fSolidMap [pSolid];
      fPODLList.append (fDisplayListId);
      fPODLTransformList.append (*fpObjectTransformation);
    }
    else {
      G4VSceneHandler::RequestPrimitives (solid);
      fSolidMap [pSolid] = fDisplayListId;
    }
  }
}

G4int G4OpenGLStoredSceneHandler::fSceneIdCount = 0;

G4int G4OpenGLStoredSceneHandler::fSceneCount = 0;

#endif
