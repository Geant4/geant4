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
// $Id: geomdefs.hh,v 1.3 2001-07-11 09:59:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// Constants, typedefs, enums for Geometry Section
//
// History:
// 30.06.95 P.Kent

#ifndef GEOMDEFS_HH
#define GEOMDEFS_HH

#include "globals.hh"

// `Infinity' - Distance returned for no intersection etc.
const G4double kInfinity = 9E99;

// Thickness of shapes for Inside function / tracking.
// Should be greater than largest math error from the shape 
// distance calculation routines.
// Tolerance is centred on surface: Inside routine uses a
//                                  tolerance dx +/- kTol/2
// Note: values not `tuned', and because of approximations kRadtolerance and
//       kAngTolerance may not always be used as an exact radius
const G4double kCarTolerance = 1E-9;
const G4double kRadTolerance = 1E-9;
const G4double kAngTolerance = 1E-9;

// Define axes for function params etc.
// X/Y/ZAxis = Normal Catesian axes
// Radial2D = Radial axis in cylindrical polar
// Radial3D = Radial axis in spherical polar
enum EAxis {kXAxis,kYAxis,kZAxis,kRadial2D,kRadial3D};

// VShape::Inside routine return codes
// kSurface => within tolerance of exact surface
enum EInside {kOutside,kSurface,kInside};

#endif





