// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateViewer.hh,v 1.5 2001-02-23 15:43:07 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLImmediateViewer : Encapsulates the `immediateness' of
//                                 an OpenGL viewer, for inheritance by
//                                 derived (X, Xm...) classes.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLIMMEDIATEVIEWER_HH
#define G4OPENGLIMMEDIATEVIEWER_HH

#include "G4VViewer.hh"
#include "G4OpenGLViewer.hh"
#include "G4OpenGLSceneHandler.hh"
#include "G4OpenGLImmediateSceneHandler.hh"
#include "G4OpenGLTransform3D.hh"
#include "globals.hh"

class G4OpenGLSceneHandler;
class G4OpenGLImmediateSceneHandler;

class G4OpenGLImmediateViewer: virtual public G4OpenGLViewer {
  
public:
  G4OpenGLImmediateViewer (G4OpenGLImmediateSceneHandler& scene);
  
private:
  G4OpenGLImmediateSceneHandler&          fSceneHandler; // Graphics Scene for this view.
};

#endif

#endif
