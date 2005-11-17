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
// $Id: testG4TwistedTrap.cc,v 1.6 2005-11-17 16:59:44 link Exp $
// GEANT4 tag $Name: 
//

// testG4TwistedTrap
//
//  Test file for class G4TwistedTrap
//
//             Ensure asserts are compiled in

#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"

#include "G4ThreeVector.hh"
#include "G4TwistedTrd.hh"
#include "G4TwistedTrap.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

const G4String OutputInside(const EInside a)
{
	switch(a) 
        {
		case kInside:  return "Inside"; 
		case kOutside: return "Outside";
		case kSurface: return "Surface";
	}
	return "????";
}


G4bool testG4TwistedTrap()
{
    G4ThreeVector pzero(0,0,0);
    G4ThreeVector pout(100*cm,100*cm,100*cm) ;
    G4ThreeVector dir1 = pzero-pout ;
    
    dir1 *= 1./dir1.mag();    

    G4TwistedTrd t1("Solid Twisted Trd #1",
		     10*cm,
		     15*cm,
		     20*cm,
		     20*cm,
		     80*cm,
		     60*deg);


    G4TwistedTrap t2("Solid Twisted Regular Trap",
		     60*deg,
		     10*cm,
		     15*cm,
		     20*cm,
		     80*cm) ;

    G4TwistedTrap t3("Solid General Twisted Trap",
		     60*deg,
		     80*cm,
		     10*deg,
		     40*deg,
		     16*cm,
		     24*cm,
		     14*cm,
		     8*cm,
		     16*cm,
		     11*cm,
		     50*deg) ;

    G4double dist = t1.DistanceToIn(pout,dir1) ;
    G4cout << "distance to = " << t1.GetName() << " = " << dist/cm << " cm" << G4endl ;

    dist = t2.DistanceToIn(pout,dir1) ;
    G4cout << "distance to = " << t2.GetName() << " = " << dist/cm << " cm" << G4endl ;

    dist = t3.DistanceToIn(pout,dir1) ;
    G4cout << "distance to = " << t3.GetName() << " = " << dist/cm << " cm" << G4endl ;


// Check Inside

    G4cout << "Test " << t1.GetName() << " Inside(pzero) is "   << t1.Inside(pzero) << G4endl ;
    G4cout << "Test " << t1.GetName() << " Inside(pout) is "    << t1.Inside(pout) << G4endl ;
    G4cout << "Test " << t2.GetName() << " Inside(pzero) is "   << t2.Inside(pzero) << G4endl ;
    G4cout << "Test " << t2.GetName() << " Inside(pout) is "    << t2.Inside(pout) << G4endl ;

    G4cout << "Test " << t3.GetName() << " Inside(pzero) is "   << t3.Inside(pzero) << G4endl ;
    G4cout << "Test " << t3.GetName() << " Inside(pout) is "    << t3.Inside(pout) << G4endl ;
    //assert(t1.Inside(pzero)==kInside) ;

    // test GetPointOnSurface

    G4ThreeVector Spoint ;
    G4ThreeVector Opoint ;
    G4ThreeVector dir ;

    EInside side ;
    for ( int i = 0 ; i < 10 ; i++ ) {
      //  G4cout << "Event " << i << G4endl << G4endl ;
      Spoint = t3.GetPointOnSurface() ;
      Opoint = t3.GetPointInSolid(Spoint.z()) ;
      dir = Opoint - Spoint ;
      
      side = t3.Inside(Spoint) ;
      //      G4cout << "SInside = " << OutputInside(side) << G4endl ;
      dist = t3.DistanceToIn(Spoint,dir/dir.mag()) ;
      G4cout << "Spoint " << Spoint << " " <<  dist << G4endl ;
    }

    return true;
}

int main()
{
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4TwistedTrap());
    return 0;
}


