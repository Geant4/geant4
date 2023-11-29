//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
#ifndef G4TRANSFORM3D_HH
#define G4TRANSFORM3D_HH

#include <CLHEP/Geometry/Transform3D.h>

using G4Transform3D = HepGeom::Transform3D;

using G4Rotate3D = HepGeom::Rotate3D;
using G4RotateX3D = HepGeom::RotateX3D;
using G4RotateY3D = HepGeom::RotateY3D;
using G4RotateZ3D = HepGeom::RotateZ3D;

using G4Translate3D = HepGeom::Translate3D;
using G4TranslateX3D = HepGeom::TranslateX3D;
using G4TranslateY3D = HepGeom::TranslateY3D;
using G4TranslateZ3D = HepGeom::TranslateZ3D;

using G4Reflect3D = HepGeom::Reflect3D;
using G4ReflectX3D = HepGeom::ReflectX3D;
using G4ReflectY3D = HepGeom::ReflectY3D;
using G4ReflectZ3D = HepGeom::ReflectZ3D;

using G4Scale3D = HepGeom::Scale3D;
using G4ScaleX3D = HepGeom::ScaleX3D;
using G4ScaleY3D = HepGeom::ScaleY3D;
using G4ScaleZ3D = HepGeom::ScaleZ3D;

#endif /* G4TRANSFORM3D_HH */
