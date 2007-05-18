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
// $Id: testSolidComparisons.cc,v 1.6 2007-05-18 11:06:34 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
// 
// Test for comparison of solids of different type but similar topology
// Returns 0 if there no inconsistencies in the answers provided by the
// compared solids.
//
// Author: Dionysios Anninos
//
// --------------------------------------------------------------
#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"
#include "G4GeometryTolerance.hh"

#include "G4ThreeVector.hh"
#include "G4TwoVector.hh"

#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4EllipticalTube.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Box.hh"
#include "G4Trap.hh"
#include "G4ExtrudedSolid.hh"

#include "Randomize.hh"

#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

using namespace CLHEP;

void logErrors(G4double x, G4double y, G4double z, 
		 G4double vx, G4double vy, G4double vz, G4double dist)
{
  G4cout <<"The point ("<<x<<","<<y<<","<<z<<")"<<G4endl
	 <<"In the direction ("<<vx<<","<<vy<<","<<vz<<")"<<G4endl
         <<"Gives a difference of "<<dist<<"cm !"<<G4endl;
}

G4bool compareEllipsoidtoOrb(G4int N)
{
  
  G4bool what = true;  
  G4int i=0, n=0;
  G4ThreeVector pin, pout, dir;
  G4double xin, yin, zin, xout, yout, zout, dist1, dist2, dist;

  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  // construct the ellipsoid and Orb with the same dimensions 
  
  G4Ellipsoid t1("Solid Ellipsoid #1",
		 20*cm,       // xSemiAxis
		 20*cm,       // ySemiAxis
		 20*cm) ;     // zSemiAxis
  
  
  G4Orb t2("Solid Orb #1", 20*cm) ; 
 
  for(i=0; i<N; i++)
  {
    xout = RandFlat::shoot(25.0*cm,300.0*cm);
    yout = RandFlat::shoot(25.0*cm,300.0*cm);
    zout = RandFlat::shoot(25.0*cm,300.0*cm);
    
    xin  = RandFlat::shoot(-10.0*cm,10.0*cm);
    yin  = RandFlat::shoot(-10.0*cm,10.0*cm);
    zin  = RandFlat::shoot(-10.0*cm,10.0*cm);
    
    pin  = G4ThreeVector(xin,  yin,  zin );
    pout = G4ThreeVector(xout, yout, zout);
    
    dir  = pin - pout;
    dir /= dir.mag();

    dist1 = t1.DistanceToIn(pout,dir) ;
    dist2 = t2.DistanceToIn(pout,dir) ;
 
    if(dist1 != kInfinity && dist2 !=kInfinity)
    {
      if(std::fabs(dist1 - dist2) >= 5.*kCarTolerance)
      { 
	what = false; 
	dist = std::fabs(dist1 - dist2);
	logErrors( pout.x(),   pout.y(),   pout.z(),
		   dir.x() ,   dir.y() ,   dir.z() , dist);
	n++;			
      } 
    }

    dist1 = t1.DistanceToOut(pin,dir);
    dist2 = t2.DistanceToOut(pin,dir);
 
    if(dist1 != kInfinity && dist2 !=kInfinity)
    {
      if(std::fabs(dist1 - dist2) >= 5.*kCarTolerance)
      {
	what = false; 
	dist = std::fabs(dist1 - dist2);
	logErrors( pin.x(), pin.y(), pin.z(),
		   dir.x(), dir.y(), dir.z(), dist);
	n++;	   
      }
    }
  }
  
  G4cout <<"The number of inconsistencies when comparing Ellipsoid to Orb were: "<<n<<"."<<G4endl;
  return what;
}

G4bool compareEllipsoidtoSphere(G4int N)
{
  
  G4bool what = true;  
  G4int i=0,n=0;
  G4ThreeVector pin, pout, dir;
  G4double xin, yin, zin, xout, yout, zout, dist1, dist2, dist;
    
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  // construct the ellipsoid and sphere with the same dimensions 
  
  G4Ellipsoid  t1("Solid Ellipsoid #1",
				   20*cm,       // xSemiAxis
				   20*cm,       // ySemiAxis
				   20*cm) ;     // zSemiAxis
  
  
  G4Sphere  t2("Solid Sphere #1", 
			     0*cm,
			     20*cm,
			     0*rad, 2*pi*rad,
			     0*rad, 2*pi*rad) ;  
 
  for(i=0; i<N; i++)
  {
    xout = RandFlat::shoot(25.0*cm,300.0*cm);
    yout = RandFlat::shoot(25.0*cm,300.0*cm);
    zout = RandFlat::shoot(25.0*cm,300.0*cm);
    
    xin  = RandFlat::shoot(-10.0*cm,10.0*cm);
    yin  = RandFlat::shoot(-10.0*cm,10.0*cm);
    zin  = RandFlat::shoot(-10.0*cm,10.0*cm);
    
    pin  = G4ThreeVector(xin, yin, zin);
    pout = G4ThreeVector(xout, yout, zout);
    
    dir  = pin - pout;
    dir /= dir.mag();

    dist1 = t1.DistanceToIn(pout,dir) ;
    dist2 = t2.DistanceToIn(pout,dir) ;
 
    if(dist1 != kInfinity && dist2 != kInfinity)
    {
      if(std::fabs(dist1 - dist2) >= 5.*kCarTolerance)
      { 
	what = false; 
	dist = std::fabs(dist1 - dist2);
	logErrors(pout.x(), pout.y(), pout.z(),
		  dir.x() , dir.y() , dir.z() , dist);
	n++;
      }
    }

    dist1 = t1.DistanceToOut(pin,dir);
    dist2 = t2.DistanceToOut(pin,dir);
 
    if(dist1 != kInfinity && dist2 !=kInfinity)
    {
      if(std::fabs(dist1 - dist2) >= 5.*kCarTolerance)
      {
	what = false; 
	dist = std::fabs(dist1 - dist2);
	logErrors( pin.x(), pin.y(), pin.z(),
		   dir.x(), dir.y(), dir.z(), dist);
	n++;
      }
    }
  }
  
  G4cout <<"The number of inconsistencies when comparing Ellipsoid to Sphere were: "<<n<<"."<<G4endl;
  return what;
}

G4bool compareEllipticalTubetoTubs(G4int N)
{
  
  G4bool what = true;  
  G4int i=0,n=0;
  G4ThreeVector pin, pout, dir;
  G4double xin, yin, zin, xout, yout, zout, dist1, dist2,dist;
    
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  // construct the tube and elliptical tube with the same dimensions 
  
  G4EllipticalTube t1("Solid EllipticalTube #1",
		      20*cm,       // xSemiAxis
		      20*cm,       // ySemiAxis
		      20*cm) ;     // zHeight
  
  
  G4Tubs t2("Solid Tubs #1", 
	    0*cm,
	    20*cm,
	    20*cm,
	    0., twopi); 
  
  for(i=0; i<N; i++)
  {
    xout = RandFlat::shoot(25.0*cm,300.0*cm);
    yout = RandFlat::shoot(25.0*cm,300.0*cm);
    zout = RandFlat::shoot(25.0*cm,300.0*cm);
    
    xin  = RandFlat::shoot(-19.0*cm,19.0*cm);
    yin  = RandFlat::shoot(-1.0*cm , 1.0*cm)*std::sqrt(361.*cm*cm-sqr(xin));
    zin  = RandFlat::shoot(-19.0*cm,19.0*cm);
    
    pin  = G4ThreeVector(xin, yin, zin);
    pout = G4ThreeVector(xout, yout, zout);
    
    dir  = pin - pout;
    dir /= dir.mag();

    dist1 = t1.DistanceToIn(pout,dir) ;
    dist2 = t2.DistanceToIn(pout,dir) ;
 
    if(dist1 != kInfinity && dist2 != kInfinity)
    {
      if(std::fabs(dist1 - dist2) >= 5.*kCarTolerance)
      {  
	what = false; 
	dist = std::fabs(dist1 - dist2);
	logErrors(pout.x(), pout.y(), pout.z(),
		  dir.x() , dir.y() , dir.z() , dist);
	n++;
      }
    }

//     dist1 = t1.DistanceToOut(pin,dir);
//     dist2 = t2.DistanceToOut(pin,dir);
 
//     if(dist1 != kInfinity && dist2 !=kInfinity)
//     {
//       if(std::fabs(dist1 - dist2) >= 5.*kCarTolerance)
//       {
// 	what = false; 
// 	dist = std::fabs(dist1 - dist2);
// 	logErrors(pin.x(), pin.y(), pin.z(),
// 		  dir.x(), dir.y(), dir.z(), dist);
// 	n++;
//       }
//     }
  }
  
  G4cout <<"The number of inconsistencies when comparing EllipticalTube to Tubs were: "<<n<<"."<<G4endl;
  return what;
}

G4bool compareBoxtoTrap(G4int N)
{
  
  G4bool what = true;  
  G4int i=0,n=0;
  G4ThreeVector pin, pout, dir;
  G4double xin, yin, zin, xout, yout, zout, dist1, dist2, dist;
    
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  G4ThreeVector pt[8] = { G4ThreeVector(-20*cm,-20*cm,-20*cm ),
			  G4ThreeVector( 20*cm,-20*cm,-20*cm ),
			  G4ThreeVector(-20*cm, 20*cm,-20*cm ),
			  G4ThreeVector( 20*cm, 20*cm,-20*cm ),
			  G4ThreeVector(-20*cm,-20*cm, 20*cm ),
			  G4ThreeVector( 20*cm,-20*cm, 20*cm ),
			  G4ThreeVector(-20*cm, 20*cm, 20*cm ),
			  G4ThreeVector( 20*cm, 20*cm, 20*cm )  };
  
  // construct the Trap  and Box with the same dimensions 
  
  G4Trap  t1("Solid Trap #1",
	     pt) ;     
  
  
  G4Box  t2("Solid Box #1", 
	    20*cm,
	    20*cm,
	    20*cm); 
  
  for(i=0; i<N; i++)
  {
    xout = RandFlat::shoot(25.0*cm,300.0*cm);
    yout = RandFlat::shoot(25.0*cm,300.0*cm);
    zout = RandFlat::shoot(25.0*cm,300.0*cm);
    
    xin  = RandFlat::shoot(-10.0*cm,10.0*cm);
    yin  = RandFlat::shoot(-10.0*cm,10.0*cm);
    zin  = RandFlat::shoot(-10.0*cm,10.0*cm);
    
    pin  = G4ThreeVector(xin, yin, zin);
    pout = G4ThreeVector(xout, yout, zout);
    
    dir  = pin - pout;
    dir /= dir.mag();

    dist1 = t1.DistanceToIn(pout,dir) ;
    dist2 = t2.DistanceToIn(pout,dir) ;
 
    if(dist1 != kInfinity && dist2 !=kInfinity)
    {
      if(std::fabs(dist1 - dist2) >= 5.*kCarTolerance)
      { 
	what = false; 
	dist = std::fabs(dist1 - dist2);
	logErrors(pout.x(),pout.y(),pout.z(),
		  dir.x(), dir.y(), dir.z(), dist);
	n++;
      }
    }

    dist1 = t1.DistanceToOut(pin,dir);
    dist2 = t2.DistanceToOut(pin,dir);
 
    if(dist1 != kInfinity && dist2 !=kInfinity)
    {
      if(std::fabs(dist1 - dist2) >= 5.*kCarTolerance)
      {
	what = false; 
	dist = std::fabs(dist1 - dist2);
	logErrors(pin.x(), pin.y(), pin.z(),
		  dir.x(), dir.y(), dir.z(), dist);
	n++;
      }
    }
  }
  
  G4cout <<"The number of inconsistencies when comparing Box to Trap were: "<<n<<"."<<G4endl;
  return what;
}

G4bool compareSpheretoOrb(G4int N)
{
  
  G4bool what = true;  
  G4int i=0,n=0;
  G4ThreeVector pin, pout, dir;
  G4double xin, yin, zin, xout, yout, zout, dist1, dist2,dist;
    
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  // construct the ellipsoid and sphere with the same dimensions 
  
  G4Orb  t1("Solid Orb #1",
		 20*cm) ;

  G4Sphere t2("Solid Sphere #1", 
			     0*cm,
			     20*cm,
			     0*deg, 2*pi*rad,
			     0*rad, 2*pi*rad) ; 
  
  for(i=0; i<N; i++)
  {
    xout = RandFlat::shoot(25.0*cm,300.0*cm);
    yout = RandFlat::shoot(25.0*cm,300.0*cm);
    zout = RandFlat::shoot(25.0*cm,300.0*cm);
    
    xin  = RandFlat::shoot(-10.0*cm,10.0*cm);
    yin  = RandFlat::shoot(-10.0*cm,10.0*cm);
    zin  = RandFlat::shoot(-10.0*cm,10.0*cm);
    
    pin  = G4ThreeVector(xin, yin, zin);
    pout = G4ThreeVector(xout, yout, zout);
    
    dir  = pin - pout;
    dir /= dir.mag();

    dist1 = t1.DistanceToIn(pout,dir) ;
    dist2 = t2.DistanceToIn(pout,dir) ;
  
    if(dist1 != kInfinity && dist2 != kInfinity)
    {
      if(std::fabs(dist1 - dist2) >= 5.*kCarTolerance)
      { 
	what = false; 
	dist = std::fabs(dist1 - dist2);
	logErrors(pout.x(),pout.y(), pout.z(),
		  dir.x() ,dir.y() , dir.z() ,dist);
	n++;
      }
    }

    dist1 = t1.DistanceToOut(pin,dir);
    dist2 = t2.DistanceToOut(pin,dir);
 
    if(dist1 != kInfinity && dist2 !=kInfinity)
    {
      if(std::fabs(dist1 - dist2) >= 5.*kCarTolerance)
      {
	what = false; 
	dist = std::fabs(dist1 - dist2);
	logErrors(pin.x(), pin.y(), pin.z(),
		  dir.x(), dir.y(), dir.z(), dist);
	n++;
      }
    }
  }
  G4cout <<"The number of inconsistencies when comparing Sphere to Orb were: "<<n<<"."<<G4endl;
  return what;
}

G4bool compareBoxtoExtruded(G4int N)
{  

  G4bool what = true;  
  G4int i=0,n=0;
  G4ThreeVector pin, pout, dir;
  G4double xin, yin, zin, xout, yout, zout, dist1, dist2, dist;
 
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  // construct the Extruded-Solid and Box with the same dimensions 
  
  std::vector<G4TwoVector> polygon;
  polygon.push_back(G4TwoVector(-1.0*cm,-0.5*cm));
  polygon.push_back(G4TwoVector(-1.0*cm, 0.5*cm));
  polygon.push_back(G4TwoVector( 1.0*cm, 0.5*cm));
  polygon.push_back(G4TwoVector( 1.0*cm,-0.5*cm));

  G4ExtrudedSolid t1("Solid Xtru #1", polygon, 3.0*cm,
                     G4TwoVector(), 1.0, G4TwoVector(), 1.0);

  G4Box  t2("Solid Box #2",
            1.0*cm,
            0.5*cm,
            3.0*cm); 
  
  for(i=0; i<N; i++)
  {
    xout = RandFlat::shoot(5.0*cm,20.0*cm);
    yout = RandFlat::shoot(5.0*cm,20.0*cm);
    zout = RandFlat::shoot(5.0*cm,20.0*cm);
    
    xin  = RandFlat::shoot(-0.9*cm,0.9*cm);
    yin  = RandFlat::shoot(-0.4*cm,0.4*cm);
    zin  = RandFlat::shoot(-2.9*cm,2.9*cm);
    
    pin  = G4ThreeVector(xin, yin, zin);
    pout = G4ThreeVector(xout, yout, zout);
    
    dir  = pin - pout;
    dir /= dir.mag();

    dist1 = t1.DistanceToIn(pout,dir) ;
    dist2 = t2.DistanceToIn(pout,dir) ;
 
    if(dist1 != kInfinity && dist2 !=kInfinity)
    {
      if(std::fabs(dist1 - dist2) >= 5.*kCarTolerance)
      { 
	what = false; 
	dist = std::fabs(dist1 - dist2);
	logErrors(pout.x(),pout.y(),pout.z(),
		  dir.x(), dir.y(), dir.z(), dist);
	n++;
      }
    }

    dist1 = t1.DistanceToOut(pin,dir);
    dist2 = t2.DistanceToOut(pin,dir);
 
    if(dist1 != kInfinity && dist2 !=kInfinity)
    {
      if(std::fabs(dist1 - dist2) >= 5.*kCarTolerance)
      {
	what = false; 
	dist = std::fabs(dist1 - dist2);
	logErrors(pin.x(), pin.y(), pin.z(),
		  dir.x(), dir.y(), dir.z(), dist);
	n++;
      }
    }
  }
  
  G4cout << "The number of inconsistencies when comparing Box to Extruded-Solid were: "
         << n << "." << G4endl;
  return what;
}

// other comparisons can be, trap with polyhedra, cube with polyhedra, tet with polyhedra

int main()
{

  G4bool what;

  what = compareEllipsoidtoOrb(1000);
  what = compareEllipticalTubetoTubs(1000);
  what = compareEllipsoidtoSphere(1000);
  what = compareBoxtoTrap(1000);
  what = compareSpheretoOrb(1000);
  what = compareBoxtoExtruded(1000);


  return 0;
}
