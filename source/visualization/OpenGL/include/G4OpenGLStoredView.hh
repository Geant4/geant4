// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredView.hh,v 1.1 1999-01-07 16:14:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLStoredView : Encapsulates the `storedness' of
//                            an OpenGL view, for inheritance by
//                            derived (X, Xm...) classes.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLSTOREDVIEW_HH
#define G4OPENGLSTOREDVIEW_HH

#include "G4VView.hh"
#include "G4OpenGLView.hh"
#include "G4OpenGLScene.hh"
#include "G4OpenGLStoredScene.hh"
#include "G4OpenGLTransform3D.hh"
#include "globals.hh"
#include <rw/tvordvec.h>

class G4OpenGLStoredScene;

class G4OpenGLStoredView: virtual public G4OpenGLView {
  
public:
  G4OpenGLStoredView (G4OpenGLStoredScene& scene);
  
protected:
  void KernelVisitDecision ();
  void DrawDisplayLists ();
  G4OpenGLStoredScene&            fScene; // Graphics Scene for this view.
};

#endif

#endif
