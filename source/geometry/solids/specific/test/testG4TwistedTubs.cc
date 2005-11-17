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
// $Id: testG4TwistedTubs.cc,v 1.3 2005-11-17 16:59:45 link Exp $
// GEANT4 tag $Name: 
//

// testG4TwistedTubs
//
//  Test file for class G4TwistedTubs
//
//             Ensure asserts are compiled in

#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4TwistedTubs.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

G4bool testG4TwistedTubs()
{
    G4ThreeVector pzero(0,0,0);

    G4ThreeVector pbigx(100,0,0),pbigy(0,100,0),pbigz(0,0,100);
    G4ThreeVector pbigmx(-100,0,0),pbigmy(0,-100,0),pbigmz(0,0,-100);
    G4ThreeVector pin1 ;

    double Rin=10*cm;
    double Rout=15*cm;
    double alpha=10*deg;
    double phi1=90*deg;
    double len=40*cm;

    G4TwistedTubs t1("Solid Twisted Tub #1",alpha,Rin,Rout,len,phi1);    

    G4ThreeVector Spoint ;
    G4ThreeVector Opoint ;
    G4ThreeVector dir ;
    G4double dist ;
    G4ThreeVector nvec ;

    for ( int i = 0 ; i < 10 ; i++ ) {
      //  G4cout << "Event " << i << G4endl << G4endl ;
      Spoint = t1.GetPointOnSurface() ;
      nvec = t1.SurfaceNormal(Spoint) ;
      dir = - nvec ;
      dist = t1.DistanceToIn(Spoint,dir/dir.mag()) ;
      G4cout << "Spoint " << Spoint << " " <<  dist << G4endl ;
    }


    return true;
}

int main()
{
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4TwistedTubs());
    return 0;
}

