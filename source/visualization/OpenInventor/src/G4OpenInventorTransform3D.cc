// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenInventorTransform3D.cc,v 1.1 1999-01-07 16:15:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// jck 17 Dec 1996
// G4OpenInventorTransform3D provides OpenGL style transformation matrix
// from G4Transform3D.

#ifdef G4VIS_BUILD_OI_DRIVER

#include <Inventor/fields/SoSFMatrix.h>

#include "G4OpenInventorTransform3D.hh"

G4OpenInventorTransform3D::G4OpenInventorTransform3D (const G4Transform3D &t) 
: G4Transform3D (t) {
  m[0]  = xx; 
  m[1]  = yx; 
  m[2]  = zx; 
  m[3]  = 0;
  m[4]  = xy; 
  m[5]  = yy; 
  m[6]  = zy; 
  m[7]  = 0;
  m[8]  = xz; 
  m[9]  = yz; 
  m[10] = zz; 
  m[11] = 0;
  m[12] = dx; 
  m[13] = dy; 
  m[14] = dz; 
  m[15] = 1;
}

SoSFMatrix* G4OpenInventorTransform3D::GetOIMatrix () const {
  SoSFMatrix *tm = new SoSFMatrix;
  tm->setValue(m[0],m[1],m[2],m[3],
               m[4],m[5],m[6],m[7],
               m[8],m[9],m[10],m[11],
               m[12],m[13],m[14],m[15]);
  return tm;
}

#endif
