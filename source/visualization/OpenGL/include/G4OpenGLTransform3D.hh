// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLTransform3D.hh,v 1.2 1999-01-09 16:22:51 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  24th October 1996
// G4OpenGLTransform3D provides OpenGL style transformation matrix
// from G4Transform3D.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#ifndef G4OPENGLTRANSFORM3D_HH
#define G4OPENGLTRANSFORM3D_HH

#include "G4Transform3D.hh"

class G4OpenGLTransform3D : public G4Transform3D {
public:
  G4OpenGLTransform3D (const G4Transform3D &t);
  const GLdouble* GetGLMatrix ();
private:
  GLdouble m[16];
};

#endif

#endif
