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
// $Id: testG4Ellipsoid.cc
// GEANT4 tag $Name: 
//
// testG4Ellipsoid
//
//  Test file for class G4Ellipsoid
//
//             Ensure asserts are compiled in

#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"

#include "G4ThreeVector.hh"
#include "G4Ellipsoid.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"


G4bool testG4Ellipsoid()
{

    //test full ellipsoid

    G4ThreeVector pzero(0,0,0);
    G4ThreeVector pout ( 200*cm,200*cm,200*cm) ;
    G4ThreeVector dir = pzero-pout ;
    dir *= 1/dir.mag();

    G4Ellipsoid t1("Solid Ellipsoid #1",
		   20*cm,       // xSemiAxis
		   10*cm,       // ySemiAxis
		   15*cm) ;     // zSemiAxis
		           

    G4double dist = t1.DistanceToIn(pout,dir) ;
    G4cout << "distance = " << dist << G4endl ;

    // Check Inside

    G4cout << "Test Inside(p) " << G4endl << G4endl ;

    G4cout << "pzero : " << t1.Inside(pzero) << G4endl ;
   
    //test the name
    G4cout << "The name is" << t1.GetName() << G4endl ; 

    // testing the volume

    G4double volume = t1.GetCubicVolume() ;
    G4cout << G4endl << G4endl ;
    G4cout << "Test Ellipsoid #1 Volume = " << volume / cm / cm /cm << " cm^3" << G4endl ;

 
    //test cut ellipsoid
    G4ThreeVector onTop(0*cm, 0*cm, 50*cm);
    G4ThreeVector axial = onTop - pzero;
    axial *= 1/axial.mag() ; 

    G4Ellipsoid t2("Solid Ellipsoid #2",
		   20*cm,
		   10*cm,
		   15*cm,
		   -10*cm,
		   8*cm) ;
    
    dist = t2.DistanceToIn(onTop,axial) ;
    G4cout << "distance = " << dist << G4endl ;

    // Check Inside

    G4cout << "Test Inside(p) " << G4endl << G4endl ;

    G4cout << "pzero : " << t2.Inside(pzero) << G4endl ;
    
    //test the name
    G4cout << "The name is" << t2.GetName() << G4endl ; 

    // testing the volume 

    volume = t2.GetCubicVolume() ;
    G4cout << G4endl << G4endl ;
    G4cout << "Test Ellipsoid #2 Volume = " << volume / cm / cm /cm << " cm^3" << G4endl ;

    return true;
}

int main() 
{ 
#ifdef NDEBUG 
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif 
    assert(testG4Ellipsoid()); 
    return 0; 
}

