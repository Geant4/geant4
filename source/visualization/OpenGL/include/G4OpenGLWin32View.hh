// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLWin32View.hh,v 1.1 1999-01-07 16:14:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// G4OpenGLWin32View : Class to provide WindowsNT specific
//                     functionality for OpenGL in GEANT4

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#ifndef G4OPENGLWIN32VIEW_HH
#define G4OPENGLWIN32VIEW_HH

#include "G4VView.hh"
#include "G4OpenGLScene.hh"
#include "globals.hh"
#include <rw/tvordvec.h>

//Win32 includes?

#include <GL/gl.h>
#include <GL/glx.h>
#include <GL/glu.h>

class G4OpenGLScene;

class G4OpenGLWin32View: virtual public G4OpenGLView {

public:
  G4OpenGLWin32View (G4OpenGLScene& scene);
  ~G4OpenGLWin32View ();
  void FinishView ();

protected:
  void GetWin32Connection ();
  void CreateGLWin32Context ();
  virtual void CreateMainWindow ();
  G4int WinSize_x, WinSize_y;

};

#endif

#endif
