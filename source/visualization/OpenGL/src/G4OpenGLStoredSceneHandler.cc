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
// $Id: G4OpenGLStoredSceneHandler.cc,v 1.14 2002-02-24 01:48:11 johna Exp $
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

#include "G4OpenGLStoredSceneHandler.hh"

#include "G4VPhysicalVolume.hh"

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
  if (fpViewer -> GetViewParameters ().GetDrawingStyle() == G4ViewParameters::hlr) {
    initialize_hlr = true;
  }
  ClearStore();  // ...and all that goes with it.
}

void G4OpenGLStoredSceneHandler::EndModeling () {
  // Make a List which calls the other lists.
  fTopPODL = glGenLists (1);
  if (!fTopPODL) {
    G4cout <<
      "ERROR: G4OpenGLStoredSceneHandler::EndModeling: Failure to allocate"
      "  display List for fTopPODL - suggest trying OpenGL Immediated mode."
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

    if (fpViewer -> GetViewParameters ().GetDrawingStyle() == G4ViewParameters::hlr) {
      initialize_hlr = false;
      //    glDisable (GL_POLYGON_OFFSET_FILL);
    }

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
    G4VSceneHandler::RequestPrimitives (solid);
  }
  else {
    // Stop-gap solution for display List re-use.  A proper
    // implementation would use geometry hierarchy.
    const G4VSolid* pSolid = &solid;
    if (!(fpCurrentPV -> IsReplicated ()) &&
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

G4int G4OpenGLStoredSceneHandler::fSceneCount = 0;

#endif
