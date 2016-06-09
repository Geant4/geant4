//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4Transform3D.hh,v 1.4 2005/11/04 08:18:51 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
#ifndef G4TRANSFORM3D_HH
#define G4TRANSFORM3D_HH

#include "globals.hh"
#include <CLHEP/Geometry/Transform3D.h>

typedef HepGeom::Transform3D G4Transform3D;

typedef HepGeom::Rotate3D G4Rotate3D;
typedef HepGeom::RotateX3D G4RotateX3D;
typedef HepGeom::RotateY3D G4RotateY3D;
typedef HepGeom::RotateZ3D G4RotateZ3D;

typedef HepGeom::Translate3D G4Translate3D;
typedef HepGeom::TranslateX3D G4TranslateX3D;
typedef HepGeom::TranslateY3D G4TranslateY3D;
typedef HepGeom::TranslateZ3D G4TranslateZ3D;

typedef HepGeom::Reflect3D G4Reflect3D;
typedef HepGeom::ReflectX3D G4ReflectX3D;
typedef HepGeom::ReflectY3D G4ReflectY3D;
typedef HepGeom::ReflectZ3D G4ReflectZ3D;

typedef HepGeom::Scale3D G4Scale3D;
typedef HepGeom::ScaleX3D G4ScaleX3D;
typedef HepGeom::ScaleY3D G4ScaleY3D;
typedef HepGeom::ScaleZ3D G4ScaleZ3D;

#endif /* G4TRANSFORM3D_HH */
