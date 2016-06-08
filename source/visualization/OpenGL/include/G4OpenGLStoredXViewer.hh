// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredXViewer.hh,v 1.2.8.1 1999/12/07 20:53:18 gunter Exp $
// GEANT4 tag $Name: geant4-01-00 $
//
// 
// Andrew Walkden  7th February 1997
// Class G4OpenGLStoredXViewer : a class derived from G4OpenGLXViewer and
//                               G4OpenGLStoredViewer.

#ifdef G4VIS_BUILD_OPENGLX_DRIVER

#ifndef G4OPENGLSTOREDXVIEWER_HH
#define G4OPENGLSTOREDXVIEWER_HH

#include "G4VViewer.hh"
#include "G4OpenGLStoredViewer.hh"
#include "G4OpenGLXViewer.hh"

class G4OpenGLStoredSceneHandler;

class G4OpenGLStoredXViewer:
public G4OpenGLXViewer, public G4OpenGLStoredViewer{
  
public:
  G4OpenGLStoredXViewer (G4OpenGLStoredSceneHandler& scene, const G4String& name = "");
  virtual ~G4OpenGLStoredXViewer ();
  void DrawView ();
};

#endif

#endif

