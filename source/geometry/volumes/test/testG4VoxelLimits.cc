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
// $Id$
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
    G4Exception("main","000",FatalException,"FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4VoxelLimits());
    return 0;
}

