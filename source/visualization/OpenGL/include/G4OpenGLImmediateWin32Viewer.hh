// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateWin32Viewer.hh,v 1.1 1999-01-09 16:22:37 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Class G4OpenGLImmediateWin32Viewer : a class derived from
//   G4OpenGLWin32Viewer and G4OpenGLImmediateViewer.

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#ifndef G4OpenGLIMMEDIATEWIN32VIEWER_HH
#define G4OpenGLIMMEDIATEWIN32VIEWER_HH

#include "G4VViewer.hh"
#include "G4OpenGLImmediateViewer.hh"
#include "G4OpenGLWin32Viewer.hh"

#include "globals.hh"
#include <rw/tvordvec.h>

class G4OpenGLImmediateSceneHandler;

class G4OpenGLImmediateWin32Viewer:
public G4OpenGLWin32Viewer, public G4OpenGLImmediateViewer{
  
public:
  G4OpenGLImmediateWin32Viewer (G4OpenGLImmediateSceneHandler& scene);
  void DrawView ();
};

#endif

#endif
