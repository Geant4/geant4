// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenGLTransform3D.cc,v 1.3 1999-12-15 14:54:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Andrew Walkden  24th October 1996
// G4OpenGLTransform3D provides OpenGL style transformation matrix
// from G4Transform3D.

#ifdef G4VIS_BUILD_OPENGL_DRIVER

#include <GL/gl.h>

#include "G4OpenGLTransform3D.hh"

G4OpenGLTransform3D::G4OpenGLTransform3D (const G4Transform3D &t):
  G4Transform3D (t) {}

const GLdouble* G4OpenGLTransform3D::GetGLMatrix () 
{
  m[0]  = (GLdouble)xx;
  m[1]  = (GLdouble)yx;
  m[2]  = (GLdouble)zx;
  m[3]  = (GLdouble)0;
  m[4]  = (GLdouble)xy;
  m[5]  = (GLdouble)yy;
  m[6]  = (GLdouble)zy;
  m[7]  = (GLdouble)0;
  m[8]  = (GLdouble)xz;
  m[9]  = (GLdouble)yz;
  m[10] = (GLdouble)zz;
  m[11] = (GLdouble)0;
  m[12] = (GLdouble)dx;
  m[13] = (GLdouble)dy;
  m[14] = (GLdouble)dz;
  m[15] = (GLdouble)1;

  return m;
}

#endif
