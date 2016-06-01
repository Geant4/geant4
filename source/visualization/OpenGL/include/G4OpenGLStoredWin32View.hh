// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredWin32View.hh,v 2.0 1998/07/02 16:44:32 gunter Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Class G4OpenGLStoredWin32View : a class derived from G4OpenGLWin32View and
//                                 G4OpenGLStoredView.

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#ifndef G4OPENGLSTOREDWIN32VIEW_HH
#define G4OPENGLSTOREDWIN32VIEW_HH

#include "G4VView.hh"
#include "G4OpenGLStoredView.hh"
#include "G4OpenGLWin32View.hh"

class G4OpenGLStoredScene;

class G4OpenGLStoredWin32View:
public G4OpenGLWin32View, public G4OpenGLStoredView{
  
public:
  G4OpenGLStoredWin32View (G4OpenGLStoredScene& scene);
  void DrawView ();
};

#endif

#endif

