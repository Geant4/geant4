// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Plane.hh,v 1.3 2000-08-28 08:57:47 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4Plane
//
// Class description:
//
// Utility class, placeholder for plane's equation coefficients.

// Authors: J.Sulkimo, P.Urban.
// ----------------------------------------------------------------------
#ifndef __G4Plane
#define __G4Plane

#include "globals.hh"

class G4Plane 
{

public:  // without description

    G4Plane() : a(0), b(0), c(0), d(0) {}

public:

    G4double a,b,c,d;
};

#endif
