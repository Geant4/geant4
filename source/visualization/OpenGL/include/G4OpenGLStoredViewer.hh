// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredViewer.hh,v 1.4 1999-11-11 15:38:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLStoredViewer : Encapsulates the `storedness' of
//                              an OpenGL viewer, for inheritance by
//                              derived (X, Xm...) classes.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLSTOREDVIEWER_HH
#define G4OPENGLSTOREDVIEWER_HH

#include "G4VViewer.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLStoredSceneHandler.hh"
#include "G4OpenGLTransform3D.hh"
#include "globals.hh"
#include "g4rw/tvordvec.h"

class G4OpenGLStoredSceneHandler;

class G4OpenGLStoredViewer: virtual public G4OpenGLViewer {
  
public:
  G4OpenGLStoredViewer (G4OpenGLStoredSceneHandler& scene);
  virtual ~G4OpenGLStoredViewer ();
  
protected:
  void KernelVisitDecision ();
  void DrawDisplayLists ();
  G4OpenGLStoredSceneHandler&            fSceneHandler; // Graphics Scene for this view.
};

#endif

#endif
