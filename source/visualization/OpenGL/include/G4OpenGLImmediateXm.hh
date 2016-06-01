// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLImmediateXm.hh,v 2.1 1998/11/06 13:42:16 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 
// Andrew Walkden  10th February 1997
// OpenGL graphics system factory.

#if defined (G4VIS_BUILD_OPENGLXM_DRIVER) || defined (G4VIS_USE_OPENGLXM)

#ifndef G4OPENGLIMMEDIATEXM_HH
#define G4OPENGLIMMEDIATEXM_HH

#include "G4VGraphicsSystem.hh"

class G4OpenGLImmediateXm: public G4VGraphicsSystem {
public:
  G4OpenGLImmediateXm ();
  G4VScene* CreateScene (const G4String& name = "");
  G4VView*  CreateView  (G4VScene&, const G4String& name = "");
};

#endif

#endif
