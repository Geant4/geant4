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
// $Id: testG4TwistedTubs.cc,v 1.1 2004-05-19 15:29:18 link Exp $
// GEANT4 tag $Name: 
//

// testG4TwistedTubs
//
//  Test file for class G4TwistedTubs
//
//             Ensure asserts are compiled in

#include <assert.h>
#include <math.h>

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

    double Rin=10;
    double Rout=15;
    double alpha=10;
    double phi1=90;
    double len=40;

    G4TwistedTubs t1("Solid Twisted Tub #1",alpha,Rin,Rout,len,phi1);    

// Check name
    assert(t1.GetName()=="Solid Twisted Tub #1");

// Check Inside
    assert(t1.Inside(pzero)==kOutside);
    assert(t1.Inside(pbigx)==kOutside);
    assert(t1.Inside(pbigy)==kOutside);
    assert(t1.Inside(pbigz)==kOutside);
    assert(t1.Inside(pbigmx)==kOutside);
    assert(t1.Inside(pbigmy)==kOutside);
    assert(t1.Inside(pbigmz)==kOutside);

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

