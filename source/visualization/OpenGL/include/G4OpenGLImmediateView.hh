// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateView.hh,v 1.1 1999-01-07 16:14:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLImmediateView : Encapsulates the `immediateness' of
//                               an OpenGL view, for inheritance by
//                               derived (X, Xm...) classes.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLIMMEDIATEVIEW_HH
#define G4OPENGLIMMEDIATEVIEW_HH

#include "G4VView.hh"
#include "G4OpenGLView.hh"
#include "G4OpenGLScene.hh"
#include "G4OpenGLImmediateScene.hh"
#include "G4OpenGLTransform3D.hh"
#include "globals.hh"
#include <rw/tvordvec.h>

class G4OpenGLScene;
class G4OpenGLImmediateScene;

class G4OpenGLImmediateView: virtual public G4OpenGLView {
  
public:
  G4OpenGLImmediateView (G4OpenGLImmediateScene& scene);
  
private:
  G4OpenGLImmediateScene&          fScene; // Graphics Scene for this view.
};

#endif

#endif
