// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateXViewer.hh,v 1.1 1999-01-09 16:22:39 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLImmediateXViewer : a class derived from G4OpenGLXViewer and
//                                  G4OpenGLImmediateViewer.

#ifdef G4VIS_BUILD_OPENGLX_DRIVER

#ifndef G4OpenGLIMMEDIATEXVIEWER_HH
#define G4OpenGLIMMEDIATEXVIEWER_HH

#include "G4VViewer.hh"
#include "G4OpenGLImmediateViewer.hh"
#include "G4OpenGLXViewer.hh"

#include "globals.hh"
#include <rw/tvordvec.h>

class G4OpenGLImmediateSceneHandler;

class G4OpenGLImmediateXViewer:
public G4OpenGLXViewer, public G4OpenGLImmediateViewer{
  
public:
  G4OpenGLImmediateXViewer (G4OpenGLImmediateSceneHandler& scene,
			  const G4String& name = "");
  void DrawView ();
};

#endif

#endif
