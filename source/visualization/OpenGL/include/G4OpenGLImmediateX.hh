// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateX.hh,v 1.4 1999-05-10 14:03:41 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  10th February 1997
// OpenGL graphics system factory.

#if defined (G4VIS_BUILD_OPENGLX_DRIVER) || defined (G4VIS_USE_OPENGLX)

#ifndef G4OPENGLIMMEDIATEX_HH
#define G4OPENGLIMMEDIATEX_HH

#include "G4VGraphicsSystem.hh"

class G4OpenGLImmediateX: public G4VGraphicsSystem {
public:
  G4OpenGLImmediateX ();
  virtual ~G4OpenGLImmediateX ();
  G4VSceneHandler* CreateSceneHandler (const G4String& name = "");
  G4VViewer*  CreateViewer  (G4VSceneHandler&, const G4String& name = "");
};

#endif

#endif
