// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4VoxelLimits.cc,v 1.2 1999-12-15 14:50:29 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// testG4VoxelLimits
//             Ensure asserts are compiled in
//
// Test file for G4Voxellimit object
//
// o Each function called
// o Complete line coverage where possible
// o Main 3 logical cases for line clipper checked

#include <assert.h>

#include "globals.hh"
#include "geomdefs.hh"
#include "G4VoxelLimits.hh"
#include "G4ThreeVector.hh"

G4bool testG4VoxelLimits()
{
    G4VoxelLimits limit;
    assert(!limit.IsXLimited());
    assert(!limit.IsYLimited());
    assert(!limit.IsZLimited());

// Check line segment unchanged when no limits
    G4ThreeVector v1(-100,-100,-100),v2(100,100,100);
    assert(limit.ClipToLimits(v1,v2));
    assert(v1==G4ThreeVector(-100,-100,-100));
    assert(v2==G4ThreeVector(100,100,100));

    limit.AddLimit(kXAxis,-10,20);
    assert(limit.IsXLimited());
    assert(!limit.IsYLimited());
    assert(!limit.IsZLimited());
    assert(limit.GetMinXExtent()==-10);
    assert(limit.GetMaxXExtent()==20);

    assert(limit.ClipToLimits(v1,v2));
    assert(v1==G4ThreeVector(-10,-10,-10));
    assert(v2==G4ThreeVector(20,20,20));

    limit.AddLimit(kYAxis,-10,20);
    assert(limit.IsXLimited());
    assert(limit.IsYLimited());
    assert(!limit.IsZLimited());
    assert(limit.GetMinXExtent()==-10);
    assert(limit.GetMaxXExtent()==20);
    assert(limit.GetMinYExtent()==-10);
    assert(limit.GetMaxYExtent()==20);
    v1=G4ThreeVector(-100,-200,0);
    v2=G4ThreeVector(100,200,0);
    assert(limit.ClipToLimits(v1,v2));
    assert(v1==G4ThreeVector(-5,-10,0));
    assert(v2==G4ThreeVector(10,20,0));

    limit.AddLimit(kZAxis,-10,20);
    assert(limit.IsXLimited());
    assert(limit.IsYLimited());
    assert(limit.IsZLimited());
    assert(limit.IsLimited(kXAxis));
    assert(limit.IsLimited(kYAxis));
    assert(limit.IsLimited(kZAxis));
    assert(limit.GetMinXExtent()==-10);
    assert(limit.GetMaxXExtent()==20);
    assert(limit.GetMinYExtent()==-10);
    assert(limit.GetMaxYExtent()==20);
    assert(limit.GetMinZExtent()==-10);
    assert(limit.GetMaxZExtent()==20);

    assert(limit.GetMaxExtent(kXAxis)==20);
    assert(limit.GetMaxExtent(kYAxis)==20);
    assert(limit.GetMaxExtent(kZAxis)==20);
    assert(limit.GetMinExtent(kXAxis)==-10);
    assert(limit.GetMinExtent(kYAxis)==-10);
    assert(limit.GetMinExtent(kZAxis)==-10);

    v1=G4ThreeVector(-100,-200,-200);
    v2=G4ThreeVector(100,200,200);
    assert(limit.ClipToLimits(v1,v2));
    assert(v1==G4ThreeVector(-5,-10,-10));
    assert(v2==G4ThreeVector(10,20,20));

    return true;
}

int main()
{
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4VoxelLimits());
    return 0;
}

