// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateXViewer.hh,v 1.4 1999-12-15 14:54:03 gunter Exp $
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
#include "g4rw/tvordvec.h"

class G4OpenGLImmediateSceneHandler;

class G4OpenGLImmediateXViewer:
public G4OpenGLXViewer, public G4OpenGLImmediateViewer{
  
public:
  G4OpenGLImmediateXViewer (G4OpenGLImmediateSceneHandler& scene,
			  const G4String& name = "");
  virtual ~G4OpenGLImmediateXViewer ();
  void DrawView ();
};

#endif

#endif
