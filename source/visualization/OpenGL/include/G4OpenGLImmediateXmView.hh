// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateXmView.hh,v 1.1 1999-01-07 16:14:48 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// G4OpenGLXmView : Class derived from G4OpenGLXView, to provide
//                  (Motif) widget OpenGL functionality for GEANT4.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OPENGLIMMEDIATEXMVIEW_HH
#define G4OPENGLIMMEDIATEXMVIEW_HH

#include "G4VView.hh"
#include "G4OpenGLImmediateView.hh"
#include "G4OpenGLXmView.hh"

#include "globals.hh"
#include <rw/tvordvec.h>

class G4OpenGLImmediateScene;

class G4OpenGLImmediateXmView:
public G4OpenGLXmView, public G4OpenGLImmediateView{
  
public:
  G4OpenGLImmediateXmView (G4OpenGLImmediateScene& scene,
			   const G4String& name = "");
  void DrawView ();

};

#endif

#endif
