// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLWin32Viewer.hh,v 1.1 1999-01-09 16:22:53 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// G4OpenGLWin32Viewer : Class to provide WindowsNT specific
//                       functionality for OpenGL in GEANT4

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#ifndef G4OPENGLWIN32VIEWER_HH
#define G4OPENGLWIN32VIEWER_HH

#include "G4VViewer.hh"
#include "G4OpenGLSceneHandler.hh"
#include "globals.hh"
#include <rw/tvordvec.h>

//Win32 includes?

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

class G4OpenGLSceneHandler;

class G4OpenGLWin32Viewer: virtual public G4OpenGLViewer {

public:
  G4OpenGLWin32Viewer (G4OpenGLSceneHandler& scene);
  ~G4OpenGLWin32Viewer ();
  void FinishView ();

protected:
  void GetWin32Connection ();
  void CreateGLWin32Context ();
  virtual void CreateMainWindow ();
  G4int WinSize_x, WinSize_y;

};

#endif

#endif
