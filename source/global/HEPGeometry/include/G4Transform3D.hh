// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Transform3D.hh,v 1.2 1999-11-16 17:29:14 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4TRANSFORM3D_HH
#define G4TRANSFORM3D_HH

#include "globals.hh"
#include <CLHEP/Geometry/Transform3D.h>

typedef HepTransform3D G4Transform3D;

typedef HepRotate3D G4Rotate3D;
typedef HepRotateX3D G4RotateX3D;
typedef HepRotateY3D G4RotateY3D;
typedef HepRotateZ3D G4RotateZ3D;

typedef HepTranslate3D G4Translate3D;
typedef HepTranslateX3D G4TranslateX3D;
typedef HepTranslateY3D G4TranslateY3D;
typedef HepTranslateZ3D G4TranslateZ3D;

typedef HepReflect3D G4Reflect3D;
typedef HepReflectX3D G4ReflectX3D;
typedef HepReflectY3D G4ReflectY3D;
typedef HepReflectZ3D G4ReflectZ3D;

typedef HepScale3D G4Scale3D;
typedef HepScaleX3D G4ScaleX3D;
typedef HepScaleY3D G4ScaleY3D;
typedef HepScaleZ3D G4ScaleZ3D;

#endif /* G4TRANSFORM3D_HH */
