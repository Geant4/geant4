// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredXmViewer.hh,v 1.1 1999-01-09 16:22:50 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// Class G4OpenGLStoredXmViewer : a class derived from G4OpenGLXmViewer 
//                                and G4OpenGLStoredViewer.

#ifdef G4VIS_BUILD_OPENGLXM_DRIVER

#ifndef G4OpenGLSTOREDXMVIEWER_HH
#define G4OpenGLSTOREDXMVIEWER_HH

#include "G4VViewer.hh"
#include "G4OpenGLStoredViewer.hh"
#include "G4OpenGLXmViewer.hh"

class G4OpenGLStoredSceneHandler;

class G4OpenGLStoredXmViewer:
public G4OpenGLXmViewer, public G4OpenGLStoredViewer{
  
public:
  G4OpenGLStoredXmViewer (G4OpenGLStoredSceneHandler& scene, const G4String& name = "");
  void DrawView ();

};

#endif

#endif
