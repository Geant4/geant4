// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateXView.hh,v 2.1 1998/11/06 13:42:15 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLImmediateXView : a class derived from G4OpenGLXView and
//                                G4OpenGLImmediateView.

#ifdef G4VIS_BUILD_OPENGLX_DRIVER

#ifndef G4OpenGLIMMEDIATEXVIEW_HH
#define G4OpenGLIMMEDIATEXVIEW_HH

#include "G4VView.hh"
#include "G4OpenGLImmediateView.hh"
#include "G4OpenGLXView.hh"

#include "globals.hh"
#include <rw/tvordvec.h>

class G4OpenGLImmediateScene;

class G4OpenGLImmediateXView:
public G4OpenGLXView, public G4OpenGLImmediateView{
  
public:
  G4OpenGLImmediateXView (G4OpenGLImmediateScene& scene,
			  const G4String& name = "");
  void DrawView ();
};

#endif

#endif
