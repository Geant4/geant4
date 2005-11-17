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
// $Id: testG4Tet.cc,v 1.2 2005-11-17 14:54:10 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// testG4Tet
//
//  Test file for class G4Tet [NOT thorough]
//
//             Ensure asserts are compiled in

#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4Tet.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

G4bool testG4Tet()
{
    G4ThreeVector pzero(0,0,0);
    G4ThreeVector pnt1(10.,0.,0.),pnt2(5.0,10.,0), pnt3(5.,5.,10.);

    G4ThreeVector *pNorm,norm;
    G4bool *pgoodNorm,goodNorm;

    pNorm=&norm;
    pgoodNorm=&goodNorm;

    G4bool  goodTet;
    G4Tet   t1( "Solid Tet #1", pzero, pnt1, pnt2, pnt3, &goodTet); 

// Check name
    assert(t1.GetName()=="Solid Tet #1");

    G4ThreeVector pntA( 1.0 , 1.0 , 1.0 ); 
    G4ThreeVector pntB( 1.5 , 0.5 , 1.0 );  
    // G4ThreeVector pntBbis= 0.1 * (pnt1-pzero) + 0.1 * (pnt3-pzero); 
    G4ThreeVector pntBr023= (1.0/3.0) * (pzero + pnt2 + pnt3); 
    G4ThreeVector pntC( 0.0,  5.0 , 1.5 );  

// Check Inside
    assert(t1.Inside(pntA)==kInside);
    assert(t1.Inside(pntB)==kSurface);
    assert(t1.Inside(pntBr023)==kSurface);
    assert(t1.Inside(pntC)==kOutside);

// Check Surface Normal
    G4ThreeVector normal;
    G4ThreeVector pntOnBotSurf012( 5.0, 5.0, 0.0); 
    G4ThreeVector vmz(0,0,-1.0); 

    normal=t1.SurfaceNormal(pntOnBotSurf012);
    assert(ApproxEqual(normal,vmz));


    return true;   
}

int main()
{
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4Tet());
    return 0;
}

