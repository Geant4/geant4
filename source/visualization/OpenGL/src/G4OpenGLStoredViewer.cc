// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredViewer.cc,v 1.1 1999-01-09 16:23:19 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLStoredViewer : Encapsulates the `storedness' of
//                            an OpenGL view, for inheritance by
//                            derived (X, Xm...) classes.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#include "G4OpenGLStoredViewer.hh"

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4ios.hh"
#include <assert.h>
#include <unistd.h>

class G4OpenGLStoredSceneHandler;

G4OpenGLStoredViewer::G4OpenGLStoredViewer (G4OpenGLStoredSceneHandler& scene):
G4VViewer (scene, -1),  
G4OpenGLViewer (scene), 
fScene (scene)
{}

void G4OpenGLStoredViewer::KernelVisitDecision () {
  
  // Trigger a display List refresh if necessary.  This is a checklist
  // of relevant view parameters.

  static G4ViewParameters lastVP;  // Initialised to default.
  G4bool need = false;

  if (
      (lastVP.GetDrawingStyle ()    != fVP.GetDrawingStyle ())    ||
      (lastVP.GetRepStyle ()        != fVP.GetRepStyle ())        ||
      (lastVP.IsCulling ()          != fVP.IsCulling ())          ||
      (lastVP.IsCullingInvisible () != fVP.IsCullingInvisible ()) ||
      (lastVP.IsDensityCulling ()   != fVP.IsDensityCulling ())   ||
      (lastVP.IsCullingCovered ()   != fVP.IsCullingCovered ())   ||
      (lastVP.IsSection ()          != fVP.IsSection ())          ||
      // No need if section plane changes.
      (lastVP.IsCutaway ()          != fVP.IsCutaway ())          ||
      (lastVP.GetCutawayPlanes ().entries () !=
                              fVP.GetCutawayPlanes ().entries ()) ||
      // No need if cutaway planes change.
      (lastVP.IsExplode ()          != fVP.IsExplode ())          ||
      (lastVP.GetNoOfSides ()       != fVP.GetNoOfSides ())
      ) {
    need = true;
  }

  if (!need && lastVP.IsDensityCulling () &&
      (lastVP.GetVisibleDensity () != fVP.GetVisibleDensity ()))
    need = true;

  if (!need && lastVP.IsExplode () &&
      (lastVP.GetExplodeFactor () != fVP.GetExplodeFactor ()))
    need = true;
      
  if (need) {
    lastVP = fVP;
    NeedKernelVisit ();
  }
}

void G4OpenGLStoredViewer::DrawDisplayLists () {
  
  if (fScene.fTopPODL) glCallList (fScene.fTopPODL);
  
  G4int nTODLs = fScene.fTODLList.entries ();
  for (int i = 0; i < nTODLs; i++) {
    glPushMatrix();
    G4OpenGLTransform3D oglt (fScene.fTODLTransformList (i));
    glMultMatrixd (oglt.GetGLMatrix ());
    glCallList (fScene.fTODLList(i));
    glPopMatrix();
  }
}

#endif
