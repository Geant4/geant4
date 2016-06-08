// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PHype.cc,v 1.2 1999/11/17 10:49:02 morita Exp $
// GEANT4 tag $Name: geant4-01-01 $
//
// class G4PHype
//
// History:
// 19.06.98 A.Kimura Converted G4PHype.cc

#include "G4VSolid.hh"
#include "G4PHype.hh"
#include "G4Hype.hh"

// Constructor - check parameters, and fills protected data members
G4PHype::G4PHype(const G4Hype* theHype)
 : G4PCSGSolid(theHype->GetName())
{

    const G4double newInnerRadius = theHype->GetInnerRadius();
    const G4double newOuterRadius = theHype->GetOuterRadius();
    const G4double newInnerStereo = theHype->GetOuterStereo();
    const G4double newOuterStereo = theHype->GetInnerStereo();
    const G4double newHalfLenZ    = theHype->GetZHalfLength();

// Check z-len
    if (newHalfLenZ>0) { 
	halfLenZ=newHalfLenZ;   
    } else {
	G4Exception("Error in G4PHype::G4PHype - invalid z half-length");
    }

// Check radii
    if (newInnerRadius>=0 && newOuterRadius>=0) {
	if (newInnerRadius < newOuterRadius) { 
	    innerRadius=newInnerRadius;
	    outerRadius=newOuterRadius;
	} else { // swapping radii  (:-)
	    innerRadius=newOuterRadius;
	    outerRadius=newInnerRadius;
	}
    } else {
	G4Exception("Error in G4PHype::G4PHype - invalid radii");
    }

    innerStereo=newInnerStereo;
    outerStereo=newOuterStereo;

    // init of precalculated quantities
    tanInnerStereo2 = tan(innerStereo)*tan(innerStereo); 
    tanOuterStereo2 = tan(outerStereo)*tan(outerStereo);
    innerRadius2 = innerRadius*innerRadius;
    outerRadius2 = outerRadius*outerRadius;
    endInnerRadius2 = HypeInnerRadius2(halfLenZ);
    endOuterRadius2 = HypeOuterRadius2(halfLenZ);
    endInnerRadius = sqrt(endInnerRadius2);
    endOuterRadius = sqrt(endOuterRadius2);
}

// Destructor
G4PHype::~G4PHype()
{;}

// make a transient object
G4VSolid* G4PHype::MakeTransientObject() const {
    G4VSolid* transientObject = new G4Hype(GetName(),
					  innerRadius,
					  outerRadius,
					  innerStereo,
					  outerStereo,
					  halfLenZ);
    return transientObject;
}
