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
// $Id: testG4TwistedTrapezoid.cc,v 1.1 2004-07-29 15:13:10 link Exp $
// GEANT4 tag $Name: 
//

// testG4TwistedTrapezoid
//
//  Test file for class G4TwistedTrapezoid
//
//             Ensure asserts are compiled in

#include <assert.h>
#include <math.h>

#include "globals.hh"
#include "geomdefs.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4TwistedTrapezoid.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

G4ThreeVector Pos(G4double phi, G4double psi, G4double b, G4double L, G4double dphi ){

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


G4bool testG4TwistedTrapezoid()
{
    G4ThreeVector pzero(0,0,0);

    G4double a=40*cm;
    G4double b=80*cm;
    G4double L=160*cm;
    G4double dphi = 30.*deg;

    G4ThreeVector pin1 = Pos( dphi/2., atan(a/b), b, L, dphi ) ;
    G4cout << "pin 1 = " << pin1 << G4endl ;

    
     G4TwistedTrapezoid t1("Solid Twisted Trapezoid #1",
					dphi ,    // twisted angle
					20*cm,       // endcap inner radius
					50*cm,       // endcap outer radius
					L/2,      // half z-length
					a/2,      // half x-length
					b/2,      // half y-length
					100.*deg ); 
    

// Check name
// assert(t1.GetName()=="Solid Twisted Trapezoid #1");

// Check Inside

    // G4cout << "pzero : " << t1.Inside(pzero) << G4endl ;
    //assert(t1.Inside(pzero)==kInside) ;

    //    G4cout << "pin1 : " << t1.Inside(pin1) << G4endl ;

    //assert(t1.Inside(pin1)==kOutside);

    G4cout << "ende" << G4endl ;

    return true;
}

int main()
{
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4TwistedTrapezoid());
    return 0;
}

