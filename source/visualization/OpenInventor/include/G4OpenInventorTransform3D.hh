// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4OpenInventorTransform3D.hh,v 1.3 1999-12-15 14:54:11 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// jck  17 Dec 1996
// G4OpenGLInventorTransform3D provides SoSFMatrix transformation matrix
// from G4Transform3D.

#ifndef G4OpenInventorTRANSFORM3D_HH
#define G4OpenInventorTRANSFORM3D_HH

#ifdef G4VIS_BUILD_OI_DRIVER

#include "G4Transform3D.hh"

class SoSFMatrix;

class G4OpenInventorTransform3D : public G4Transform3D {
public:
  G4OpenInventorTransform3D (const G4Transform3D &t);
  SoSFMatrix* GetOIMatrix () const;

private:
  G4float m[16];
};

#endif

#endif
