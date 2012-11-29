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
// GEANT4 tag $Name: 
//

// testG4TwistedTubs
//
//  Test file for class G4TwistedTubs
//
// 08.10.12 V. Grichine, bug report no. 899 was added 
//
//             Ensure asserts are compiled in

#undef NDEBUG
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

///////////////////////////////////////////////////////////////////
//
// Dave's auxiliary function

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

///////////////////////////////////////////////////////////////////
//
//


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
  G4double dist;
  G4bool *pgoodNorm=NULL,calcNorm=true;
  G4ThreeVector *pNorm=NULL,norm;


  // b899 distance to out id infinity

  G4double innerRadiusOfTheTube = 407.195*mm;
  G4double outerRadiusOfTheTube = 422.92898418921311*mm;
  G4double hightOfTheTube = 1090*mm;
  //  G4double startAngleOfTheTube = 0.*deg;
  G4double spanningAngleOfTheTube = 2.25*deg;
  G4double twistedAngleOfTheTube = 0.23561944901923448*rad;

  G4TwistedTubs* tracker_tube = new G4TwistedTubs("tracker_tube",
                                                  twistedAngleOfTheTube,
                                                  innerRadiusOfTheTube,
                                                  outerRadiusOfTheTube-2.*micrometer,
                                                  hightOfTheTube,
                                                  spanningAngleOfTheTube-0.01*deg);


  G4ThreeVector PrimaryPosition;
  G4double x = -14.37245571118363*mm;
  G4double y  = -421.5761625068428*mm;
  PrimaryPosition[2] = -859.6092789812545*mm;


  spanningAngleOfTheTube = 122.*2.25*deg;

  G4double cosTheta = std::cos(spanningAngleOfTheTube);
  G4double sinTheta = std::sin(spanningAngleOfTheTube);

  // PrimaryPosition[0] = -14.37245571118363*mm;
  // PrimaryPosition[1] = -421.5761625068428*mm;
  PrimaryPosition[0] = x*cosTheta + y*sinTheta;
  PrimaryPosition[1] = -x*sinTheta + y*cosTheta;

  // G4ParticleMomentum aMomentum;
  G4ThreeVector aMomentum;

  // aMomentum[0] = 0.05348375305152415*MeV;
  // aMomentum[1] = -0.005546019422201249*MeV;
  aMomentum[2] = 0.04481714826766207*MeV;


  x = 0.05348375305152415*MeV;
  y = -0.005546019422201249*MeV;

  aMomentum[0] = x*cosTheta + y*sinTheta;
  aMomentum[1] = -x*sinTheta + y*cosTheta;

  aMomentum /= aMomentum.mag();

  // PrimaryPosition is on surface and going inside

  EInside a1 = tracker_tube->Inside(PrimaryPosition);	
  G4cout << "bug899, starting point is " << OutputInside(a1) << G4endl;
  G4ThreeVector normal;
 
  normal = tracker_tube->SurfaceNormal(PrimaryPosition);

  G4double scalar;

  scalar = aMomentum.x()*normal.x() + aMomentum.y()*normal.y() + aMomentum.z()*normal.z();

  G4cout<<"dirction*normal = "<<scalar<<G4endl; 

  dist = tracker_tube->DistanceToIn(PrimaryPosition,aMomentum);

  G4cout<<"distance to in = "<<dist<<G4endl; 


  dist = tracker_tube->DistanceToOut(PrimaryPosition,aMomentum,calcNorm,pgoodNorm,pNorm);

  G4cout<<"distance to out = "<<dist<<G4endl; 




#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4TwistedTubs());
    return 0;
}

