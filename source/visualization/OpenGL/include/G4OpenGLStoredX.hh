// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLStoredX.hh,v 1.5 1999-12-15 14:54:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// OpenGL graphics system factory.

#if defined (G4VIS_BUILD_OPENGLX_DRIVER) || defined (G4VIS_USE_OPENGLX)

#ifndef G4OPENGLSTOREDX_HH
#define G4OPENGLSTOREDX_HH

#include "G4VGraphicsSystem.hh"

class G4OpenGLStoredX: public G4VGraphicsSystem {
public:
  G4OpenGLStoredX ();
  virtual ~G4OpenGLStoredX ();
  G4VSceneHandler* CreateSceneHandler (const G4String& name = "");
  G4VViewer*  CreateViewer  (G4VSceneHandler&, const G4String& name = "");
};

#endif

#endif
