// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateWin32View.hh,v 1.1 1999-01-07 16:14:48 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Class G4OpenGLImmediateWin32View : a class derived from G4OpenGLWin32View and
//                                    G4OpenGLImmediateView.

#ifdef G4VIS_BUILD_OPENGLWIN32_DRIVER

#ifndef G4OpenGLIMMEDIATEWIN32VIEW_HH
#define G4OpenGLIMMEDIATEWIN32VIEW_HH

#include "G4VView.hh"
#include "G4OpenGLImmediateView.hh"
#include "G4OpenGLWin32View.hh"

#include "globals.hh"
#include <rw/tvordvec.h>

class G4OpenGLImmediateScene;

class G4OpenGLImmediateWin32View:
public G4OpenGLWin32View, public G4OpenGLImmediateView{
  
public:
  G4OpenGLImmediateWin32View (G4OpenGLImmediateScene& scene);
  void DrawView ();
};

#endif

#endif
