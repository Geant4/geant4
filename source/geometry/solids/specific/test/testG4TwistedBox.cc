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
// $Id: testG4TwistedBox.cc,v 1.2 2004-11-12 14:32:57 link Exp $
// GEANT4 tag $Name: 
//

// testG4TwistedBox
//
//  Test file for class G4TwistedBox
//
//             Ensure asserts are compiled in

#include <assert.h>
#include <math.h>

#include "globals.hh"
#include "geomdefs.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4TwistedBox.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

G4ThreeVector Pos(G4double phi, G4double psi, G4double b, G4double L, G4double dphi ){
  // returns a point on the surface
  G4cout << "phi = " << phi << G4endl ;
  G4cout << "psi = " << psi << G4endl ;
  G4cout << "b = " << b << G4endl ;
  G4cout << "L = " << L << G4endl ;
  G4cout << "dphi = " << dphi << G4endl ;

  G4double fx = b/2 * cos(phi) - b/2 * sin(phi) * tan (psi) ;
  G4double fy = b/2 * sin(phi) + b/2 * cos(phi) * tan (psi) ;
  G4double fz = L * phi / dphi ;
  
  G4ThreeVector vec(fx,fy,fz) ;
  return vec ;
}


G4ThreeVector NormAng( G4double phi, G4double psi, G4double b, G4double L, G4double dphi ) 
{
  // function to calculate the norm at a given point on the surface

  G4ThreeVector nvec( L*cos(phi), L*sin(phi), b*dphi*tan(psi));
  return nvec.unit() ;
}


G4bool testG4TwistedBox()
{
    G4ThreeVector pzero(0,0,0);

    G4double a=80*cm;  // y - side 
    G4double b=40*cm;  // x - side
    G4double L=160*cm;

    G4cout << "a , b, L = " << a << " , " << b << " , " << L << G4endl ;
    G4double dphi = 30.*deg;

    G4double phi1 = 10.*deg ;
    G4double psi1 = 0.1*atan(a/b);
    G4cout <<  "phi1, psi1 = " << phi1 << " , " << psi1 << G4endl ;

    G4ThreeVector psurf1 = Pos( phi1, psi1 , b, L, dphi ) ;
    G4cout << "psurf 1 = " << psurf1 << G4endl ;
    G4ThreeVector pout ( 200*cm,200*cm,200*cm) ;
    G4ThreeVector dir = psurf1-pout ;
    G4ThreeVector psurf2 ( 0,0,L/2) ;
    G4ThreeVector psurf3 ( 0,0,-L/2) ;

    G4ThreeVector pcorner1 = Pos( dphi/2, atan(a/b), b, L, dphi ) ; 
    G4ThreeVector pcorner2 = Pos( dphi/2, -atan(a/b), b, L, dphi ) ; 
    G4ThreeVector pcorner3 = Pos( -dphi/2, atan(a/b), b, L, dphi ) ; 
    G4ThreeVector pcorner4 = Pos( -dphi/2, -atan(a/b), b, L, dphi ) ; 

    G4ThreeVector normvec = NormAng(phi1,psi1,b,L,dphi) ;
    G4ThreeVector pnormout1 = psurf1 + 22.124 * normvec ;
    G4ThreeVector pnormout2 = psurf1 + 235848.124 * normvec ;
    G4double dn1 = ( psurf1 - pnormout1 ).mag() ;
    G4double dn2 = ( psurf1 - pnormout2 ).mag() ;

    // problematic points
    G4ThreeVector pTrack0 ( 88.366, 437.158, 800 ) ;
    G4ThreeVector vTrack0 ( -0.19617, 0.5389857, -0.8191519 ) ;

    G4ThreeVector pTrack1 (-74.96379937,386.6189163,800 ) ;
    G4ThreeVector vTrack1 (-0.1961747411,0.5389857236,-0.8191519156) ;

    G4ThreeVector pTrack2 (-71.41023511,386.3731183,800) ;
    G4ThreeVector vTrack2 (-0.1961747401,0.5389857236,-0.8191519158) ;


    // track with vz = 0
    G4ThreeVector pTrack3 = psurf1 + G4ThreeVector(100,0,0) ;
    G4ThreeVector vTrack3 (-10.,0.,0.) ;
   

    G4TwistedBox t1("Solid Twisted Box #1",
			  dphi ,    // twisted angle
			  b/2,       // half x-length
			  a/2,       // half y-length
			  L/2) ;      // half z-length

    G4double dist = t1.DistanceToIn(pout,dir) ;
    G4cout << "distance = " << dist << G4endl ;

    // point on the upper endcap
    G4cout << G4endl << G4endl  ;
    G4cout << "Test upper endcap ... " << G4endl << G4endl ;

    G4ThreeVector pin2 ( 0, 0, L/2 ) ;
    G4ThreeVector pout2 ( b/10., - a/10, L*2.224 ) ;
    dir = pin2 - pout2 ;
    dist = t1.DistanceToIn(pout2,dir) ;
    G4cout << "distance = " << dist << G4endl ;

// Check Inside

    G4cout << "Test Inside(p) " << G4endl << G4endl ;

    G4cout << "pzero : " << t1.Inside(pzero) << G4endl ;
    //assert(t1.Inside(pzero)==kInside) ;

    G4cout << "psurf1 : " << t1.Inside(psurf1) << G4endl ;
    G4cout << "psurf2 : " << t1.Inside(psurf2) << G4endl ;
    G4cout << "psurf3 : " << t1.Inside(psurf3) << G4endl ;

    //assert(t1.Inside(pin1)==kOutside);

    G4cout << "pout  : " << t1.Inside(pout)  << G4endl ;
    G4cout << "pout2 : " << t1.Inside(pout2) << G4endl ;
    G4cout << "pnormout1 : " << t1.Inside(pnormout1) << G4endl ;
    G4cout << "pnormout2 : " << t1.Inside(pnormout2) << G4endl ;
    G4cout << "pcorner1 : " << t1.Inside(pcorner1) << G4endl ;
    G4cout << "pcorner2 : " << t1.Inside(pcorner2) << G4endl ;
    G4cout << "pcorner3 : " << t1.Inside(pcorner3) << G4endl ;
    G4cout << "pcorner4 : " << t1.Inside(pcorner4) << G4endl ;
    G4cout << "pTrack0 : " << t1.Inside(pTrack0) << G4endl ;
    G4cout << "pTrack1 : " << t1.Inside(pTrack1) << G4endl ;
    G4cout << "pTrack2 : " << t1.Inside(pTrack2) << G4endl ;
    
    // test DistanceToIn(p)

    G4cout << G4endl << "Test DistanceToIn(p) "  << G4endl ;

    G4cout << "Distance to Pnormout1 = " << dn1 << G4endl ;
    G4cout << "pnormout1 = " << pnormout1 << G4endl ;
    G4cout << "Distance to Pnormout2 = " << dn2 << G4endl ;
    G4cout << "pnormout2 = " << pnormout2 << G4endl ;

    dist = t1.DistanceToIn(pout) ;
    G4cout << "distanceToIn pout " << pout << " = " << dist << G4endl ;
    dist = t1.DistanceToIn(pout2) ;
    G4cout << "distanceToIn pout2 " << pout2 << " = " << dist << G4endl ;
    dist = t1.DistanceToIn(pnormout1) ;
    G4cout << "distanceToIn pnormout1 " << pnormout1 << " = " << dist << G4endl ;
    dist = t1.DistanceToIn(pnormout2) ;
    G4cout << "distanceToIn pnormout2 " << pnormout2 << " = " << dist << G4endl ;
    dist = t1.DistanceToIn(psurf1) ;
    G4cout << "distanceToIn psurf1 " << psurf1 << " = " << dist << G4endl ;
    dist = t1.DistanceToIn(psurf2) ;
    G4cout << "distanceToIn psurf2 " << psurf2 << " = " << dist << G4endl ;
    dist = t1.DistanceToIn(psurf3) ;
    G4cout << "distanceToIn psurf3 " << psurf3 << " = " << dist << G4endl ;
    dist = t1.DistanceToIn(pzero) ;
    G4cout << "distanceToIn pzero " << pzero << " = " << dist << G4endl ;
    dist = t1.DistanceToIn(pcorner1) ;
    G4cout << "distanceToIn corner1 " << pcorner1 << " = " << dist << G4endl ;
    dist = t1.DistanceToIn(pcorner2) ;
    G4cout << "distanceToIn corner2 " << pcorner2 << " = " << dist << G4endl ;
    dist = t1.DistanceToIn(pcorner3) ;
    G4cout << "distanceToIn corner3 " << pcorner3 << " = " << dist << G4endl ;
    dist = t1.DistanceToIn(pcorner4) ;
    G4cout << "distanceToIn corner4 " << pcorner4 << " = " << dist << G4endl ;
    dist = t1.DistanceToIn(pTrack0) ;
    G4cout << "distanceToIn pTrack0 " << pTrack0 << " = " << dist << G4endl ;
    dist = t1.DistanceToIn(pTrack1) ;
    G4cout << "distanceToIn pTrack1 " << pTrack1 << " = " << dist << G4endl ;
    dist = t1.DistanceToIn(pTrack2) ;
    G4cout << "distanceToIn pTrack2 " << pTrack2 << " = " << dist << G4endl ;
 
    // test DistanceToOut(p)

    G4cout << G4endl << "Test DistanceToOut(p) " << G4endl ;

    dist = t1.DistanceToOut(pout) ;
    G4cout << "distanceToOut pout " << pout << " = " << dist << G4endl ;
    dist = t1.DistanceToOut(pout2) ;
    G4cout << "distanceToOut pout2 " << pout2 << " = " << dist << G4endl ;
    dist = t1.DistanceToOut(pnormout1) ;
    G4cout << "distanceToOut pnormout1 " << pnormout1 << " = " << dist << G4endl ;
    dist = t1.DistanceToOut(pnormout2) ;
    G4cout << "distanceToOut pnormout2 " << pnormout2 << " = " << dist << G4endl ;
    dist = t1.DistanceToOut(psurf1) ;
    G4cout << "distanceToOut psurf1 " << psurf1 << " = " << dist << G4endl ;
    dist = t1.DistanceToOut(psurf2) ;
    G4cout << "distanceToOut psurf2 " << psurf2 << " = " << dist << G4endl ;
    dist = t1.DistanceToOut(psurf3) ;
    G4cout << "distanceToOut psurf3 " << psurf3 << " = " << dist << G4endl ;
    dist = t1.DistanceToOut(pzero) ;
    G4cout << "distanceToOut pzero " << pzero << " = " << dist << G4endl ;
    dist = t1.DistanceToOut(pcorner1) ;
    G4cout << "distanceToOut corner1 " << pcorner1 << " = " << dist << G4endl ;
    dist = t1.DistanceToOut(pcorner2) ;
    G4cout << "distanceToOut corner2 " << pcorner2 << " = " << dist << G4endl ;
    dist = t1.DistanceToOut(pcorner3) ;
    G4cout << "distanceToOut corner3 " << pcorner3 << " = " << dist << G4endl ;
    dist = t1.DistanceToOut(pcorner4) ;
    G4cout << "distanceToOut corner4 " << pcorner4 << " = " << dist << G4endl ;
    dist = t1.DistanceToOut(pTrack0) ;
    G4cout << "distanceToOut pTrack0 " << pTrack0 << " = " << dist << G4endl ;
    dist = t1.DistanceToOut(pTrack1) ;
    G4cout << "distanceToOut pTrack1 " << pTrack1 << " = " << dist << G4endl ;
    dist = t1.DistanceToOut(pTrack2) ;
    G4cout << "distanceToOut pTrack2 " << pTrack2 << " = " << dist << G4endl ;
 

    // Test DistanceToIn(p,v)

    G4cout << G4endl << "Test DistanceToIn(p,v) " << G4endl ;

    dist = t1.DistanceToIn(pTrack0,vTrack0) ;
    G4cout << "distanceToIn(pTrack0,vTrack0) = " << dist << G4endl ;
    dist = t1.DistanceToIn(pTrack1,vTrack1) ;
    G4cout << "distanceToIn(pTrack1,vTrack1) = " << dist << G4endl ;
    dist = t1.DistanceToIn(pTrack2,vTrack2) ;
    G4cout << "distanceToIn(pTrack2,vTrack2) = " << dist << G4endl ;
    dist = t1.DistanceToIn(pTrack3,vTrack3) ;
    G4cout << "distanceToIn(pTrack3,vTrack3) = " << dist << G4endl ;

    // Normal vector
    G4ThreeVector norm ;

    G4cout << G4endl << "Test Normal " << G4endl ;

    norm = t1.SurfaceNormal(pTrack0) ;
    G4cout << "Normal vector on pTrack0 = " << norm << G4endl ;
    norm = t1.SurfaceNormal(pTrack1) ;
    G4cout << "Normal vector on pTrack1 = " << norm << G4endl ;
    norm = t1.SurfaceNormal(pTrack2) ;
    G4cout << "Normal vector on pTrack2 = " << norm << G4endl ;
    norm = t1.SurfaceNormal(psurf1) ;
    G4cout << "Normal vector on psurf1 = " << norm << G4endl ;
    norm = t1.SurfaceNormal(psurf2) ;
    G4cout << "Normal vector on psurf2 = " << norm << G4endl ;
    norm = t1.SurfaceNormal(psurf3) ;
    G4cout << "Normal vector on psurf3 = " << norm << G4endl ;

    // testing the volume

    G4double volume = t1.GetCubicVolume() ;
    G4cout << G4endl << G4endl ;
    G4cout << "Twisted Box volume = " << volume / cm / cm /cm << " cm^3" << G4endl ;
    G4cout << "should be = " << a*b*L / cm / cm / cm << " cm^3" << G4endl ;
    G4cout << "test ended" << G4endl ;

    return true;
}

int main()
{
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4TwistedBox());
    return 0;
}


