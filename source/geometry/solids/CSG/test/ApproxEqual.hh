// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ApproxEqual.hh,v 1.2 1999-12-15 14:50:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// ApproxEqual Functions for geometry test programs
//
// History:
// 20.07.95 P.Kent Translated from old code

#ifndef APPROXEQUAL_HH
#define APPROXEQUAL_HH

#include "globals.hh"
#include "geomdefs.hh"

#include "G4ThreeVector.hh"
#include "G4AffineTransform.hh"

const G4double kApproxEqualTolerance = 1E-6;

// Return true if the double check is approximately equal to target
//
// Process:
//
// Return true is check if less than kApproxEqualTolerance from target

G4bool ApproxEqual(const G4double check,const G4double target)
{
    return (fabs(check-target)<kApproxEqualTolerance) ?true:false;
}

// Return true if the 3vector check is approximately equal to target
G4bool ApproxEqual(const G4ThreeVector& check, const G4ThreeVector& target)
{
    return (ApproxEqual(check.x(),target.x())&&
	   ApproxEqual(check.y(),target.y())&&
	    ApproxEqual(check.z(),target.z())) ? true : false;
}


G4bool ApproxEqual(const G4AffineTransform &tf1,
		   const G4AffineTransform &tf2)
{
	for (G4int i=0;i<15;i++)
		{
		if (!ApproxEqual(tf1[i],tf2[i])) return false;
		}
	return true;
}

#endif









