// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredView.cc,v 1.1 1999-01-07 16:14:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLStoredView : Encapsulates the `storedness' of
//                            an OpenGL view, for inheritance by
//                            derived (X, Xm...) classes.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#include "G4OpenGLStoredView.hh"

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

#include "G4ios.hh"
#include <assert.h>
#include <unistd.h>

class G4OpenGLStoredScene;

G4OpenGLStoredView::G4OpenGLStoredView (G4OpenGLStoredScene& scene):
G4VView (scene, -1),  
G4OpenGLView (scene), 
fScene (scene)
{}

void G4OpenGLStoredView::KernelVisitDecision () {
  
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

void G4OpenGLStoredView::DrawDisplayLists () {
  
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
