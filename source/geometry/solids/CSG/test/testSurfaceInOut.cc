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
// Test File for tracking functions on a solid surface
//
// o Basic asserts on each function +
//   awkward cases for tracking / geom algorithms
//
// o Add tests on dicovering bugs in G4Sphere.cc...
//
// History:
//
// 11.11.04 V.Grichine creation for box, tubs, cons, sphere, orb, torus
// 07.06.05 J.Apostolakis:  revised to test each solid in turn,  use argc/argv

// Notes:  (J.A. June 2005)
//  - To call it use :  testSurfaceNormal  [no_points] [no_directions_each]
//  - Old 'main' is now a test function,  new main loops over each solid type

#include "G4ios.hh"
#undef NDEBUG
#include <assert.h>
#include <cmath>
#include "globals.hh"
#include "geomdefs.hh"
#include "Randomize.hh"
#include "G4GeometryTolerance.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Hype.hh"
#include "G4Para.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"



///////////////////////////////////////////////////////////////////////////
//
//


//const G4double kApproxEqualTolerance = kCarTolerance;

// Return true if the double check is approximately equal to target
//
// Process:
//
// Return true is difference < kApproxEqualTolerance

//G4bool ApproxEqual(const G4double check,const G4double target)
//{
//    return (std::fabs(check-target)<kApproxEqualTolerance) ? true : false ;
//}

// Return true if the 3vector check is approximately equal to target
//G4bool ApproxEqual(const G4ThreeVector& check, const G4ThreeVector& target)
//{
//    return (ApproxEqual(check.x(),target.x())&&
//	   ApproxEqual(check.y(),target.y())&&
//	    ApproxEqual(check.z(),target.z()))? true : false;
//}




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


/////////////////////////////////////////////////////////////////////////////
//
// Random unit direction vector

G4ThreeVector GetRandomUnitVector()
{
  G4double cosTheta, sinTheta, phi, vx, vy, vz;

  cosTheta = -1. + 2.*G4UniformRand();
  if( cosTheta > 1.)  cosTheta = 1.;
  if( cosTheta < -1.) cosTheta = -1.;
  sinTheta = std::sqrt( 1. - cosTheta*cosTheta );
  
  phi = 2*pi*G4UniformRand();

  vx = sinTheta*std::cos(phi);
  vy = sinTheta*std::sin(phi);
  vz = cosTheta;

  return G4ThreeVector(vx,vy,vz);
}

/////////////////////////////////////////////////////////////////////////////
//
// Random  vector on box surface 

G4ThreeVector GetVectorOnBox( G4Box& box )
{
  G4double rand, a, b, c, px, py, pz;
  G4double part = 1./6.;

  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  a = box.GetXHalfLength();
  b = box.GetYHalfLength();
  c = box.GetZHalfLength();
  
  rand = G4UniformRand();

  if      ( rand < part )
  {
    px = -a - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 
    py = -b - 0.5*kCarTolerance + (2.*b + kCarTolerance)*G4UniformRand(); 
    pz = -c - 0.5*kCarTolerance + (2.*c + kCarTolerance)*G4UniformRand();  
  }
  else if ( rand < 2*part )
  {
    px =  a - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 
    py = -b - 0.5*kCarTolerance + (2.*b + kCarTolerance)*G4UniformRand(); 
    pz = -c - 0.5*kCarTolerance + (2.*c + kCarTolerance)*G4UniformRand();  
  }
  else if ( rand < 3*part )
  {
    px = -a - 0.5*kCarTolerance + (2.*a + kCarTolerance)*G4UniformRand(); 
    py = -b - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 
    pz = -c - 0.5*kCarTolerance + (2.*c + kCarTolerance)*G4UniformRand();  
  }
  else if ( rand < 4*part )
  {
    px = -a - 0.5*kCarTolerance + (2.*a + kCarTolerance)*G4UniformRand(); 
    py =  b - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 
    pz = -c - 0.5*kCarTolerance + (2.*c + kCarTolerance)*G4UniformRand();  
  }
  else if ( rand < 5*part )
  {
    px = -a - 0.5*kCarTolerance + (2.*a + kCarTolerance)*G4UniformRand(); 
    py = -b - 0.5*kCarTolerance + (2.*b + kCarTolerance)*G4UniformRand(); 
    pz = -c - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();  
  }
  else 
  {
    px = -a - 0.5*kCarTolerance + (2.*a + kCarTolerance)*G4UniformRand(); 
    py = -b - 0.5*kCarTolerance + (2.*b + kCarTolerance)*G4UniformRand(); 
    pz =  c - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();  
  }

  return G4ThreeVector(px,py,pz);
}

/////////////////////////////////////////////////////////////////////////////
//
// Random vector on orb surface

G4ThreeVector GetVectorOnOrb(G4Orb& orb)
{
  G4double cosTheta, sinTheta, phi, radius, px, py, pz;

  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  radius = orb.GetRadius();
  radius += -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 

  cosTheta = -1. + 2.*G4UniformRand();
  if( cosTheta > 1.)  cosTheta = 1.;
  if( cosTheta < -1.) cosTheta = -1.;
  sinTheta = std::sqrt( 1. - cosTheta*cosTheta );
  
  phi = 2*pi*G4UniformRand();

  px = radius*sinTheta*std::cos(phi);
  py = radius*sinTheta*std::sin(phi);
  pz = radius*cosTheta;

  return G4ThreeVector(px,py,pz);
}

/////////////////////////////////////////////////////////////////////////////
//
// Random vector on sphere surface

G4ThreeVector GetVectorOnSphere(G4Sphere& sphere)
{
  G4double cosTheta, sinTheta, phi, radius, px, py, pz;
  G4double part = 1./6.;
  G4double rand = G4UniformRand();

  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4double kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();

  G4double pRmin  = sphere.GetInnerRadius();
  G4double pRmax  = sphere.GetOuterRadius();
  G4double phi1   = sphere.GetStartPhiAngle();
  G4double phi2   = phi1 + sphere.GetDeltaPhiAngle();
  G4double theta1 = sphere.GetStartThetaAngle();
  G4double theta2 = theta1 + sphere.GetDeltaThetaAngle();

  if      ( rand < part ) // Rmax
  {
    radius   = pRmax -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 

    cosTheta = std::cos(theta2+0.5*kAngTolerance) + 
              (std::cos(theta1-0.5*kAngTolerance)-std::cos(theta2+0.5*kAngTolerance))*G4UniformRand();
    if( cosTheta > 1.)  cosTheta = 1.;
    if( cosTheta < -1.) cosTheta = -1.;
    sinTheta = std::sqrt( 1. - cosTheta*cosTheta );
  
    phi      = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }
  else if ( rand < 2*part )  // Rmin
  {
    radius   = pRmin -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 

    cosTheta = std::cos(theta2+0.5*kAngTolerance) + 
              (std::cos(theta1-0.5*kAngTolerance)-std::cos(theta2+0.5*kAngTolerance))*G4UniformRand();
    if( cosTheta > 1.)  cosTheta = 1.;
    if( cosTheta < -1.) cosTheta = -1.;
    sinTheta = std::sqrt( 1. - cosTheta*cosTheta );
  
    phi      = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }
  else if ( rand < 3*part )  // phi1
  {
    radius   = pRmin - 0.5*kCarTolerance + (pRmax-pRmin+kCarTolerance)*G4UniformRand(); 

    cosTheta = std::cos(theta2+0.5*kAngTolerance) + 
              (std::cos(theta1-0.5*kAngTolerance)-std::cos(theta2+0.5*kAngTolerance))*G4UniformRand();
    if( cosTheta > 1.)  cosTheta = 1.;
    if( cosTheta < -1.) cosTheta = -1.;
    sinTheta = std::sqrt( 1. - cosTheta*cosTheta );
  
    phi      = phi1 -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();
  }
  else if ( rand < 4*part )  // phi2
  {
    radius   = pRmin - 0.5*kCarTolerance + (pRmax-pRmin+kCarTolerance)*G4UniformRand(); 

    cosTheta = std::cos(theta2+0.5*kAngTolerance) + 
              (std::cos(theta1-0.5*kAngTolerance)-std::cos(theta2+0.5*kAngTolerance))*G4UniformRand();
    if( cosTheta > 1.)  cosTheta = 1.;
    if( cosTheta < -1.) cosTheta = -1.;
    sinTheta = std::sqrt( 1. - cosTheta*cosTheta );
  
    phi      = phi2 -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();
  }
  else if ( rand < 5*part )  // theta1
  {
    radius   = pRmin - 0.5*kCarTolerance + (pRmax-pRmin+kCarTolerance)*G4UniformRand(); 

    cosTheta = std::cos(theta1+0.5*kAngTolerance) + 
              (std::cos(theta1-0.5*kAngTolerance)-std::cos(theta1+0.5*kAngTolerance))*G4UniformRand();
    if( cosTheta > 1.)  cosTheta = 1.;
    if( cosTheta < -1.) cosTheta = -1.;
    sinTheta = std::sqrt( 1. - cosTheta*cosTheta );
  
    phi      = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }
  else // theta2
  {
    radius   = pRmin - 0.5*kCarTolerance + (pRmax-pRmin+kCarTolerance)*G4UniformRand(); 

    cosTheta = std::cos(theta2+0.5*kAngTolerance) + 
              (std::cos(theta2-0.5*kAngTolerance)-std::cos(theta2+0.5*kAngTolerance))*G4UniformRand();
    if( cosTheta > 1.)  cosTheta = 1.;
    if( cosTheta < -1.) cosTheta = -1.;
    sinTheta = std::sqrt( 1. - cosTheta*cosTheta );
  
    phi      = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }

  px = radius*sinTheta*std::cos(phi);
  py = radius*sinTheta*std::sin(phi);
  pz = radius*cosTheta;

  return G4ThreeVector(px,py,pz);
}

/////////////////////////////////////////////////////////////////////////////
//
// Random vector on tubs surface

G4ThreeVector GetVectorOnTubs(G4Tubs& tubs)
{
  G4double phi, radius, px, py, pz;
  G4double part = 1./6.;
  G4double rand = G4UniformRand();

  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4double kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();

  G4double pRmin = tubs.GetInnerRadius   ();
  G4double pRmax = tubs.GetOuterRadius   ();
  G4double tubsZ = tubs.GetZHalfLength   ();
  G4double phi1  = tubs.GetStartPhiAngle ();
  G4double phi2  = phi1 + tubs.GetDeltaPhiAngle ();


  if      ( rand < part ) // Rmax
  {
    radius = pRmax -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 
    pz     = -tubsZ - 0.5*kCarTolerance + (2*tubsZ + kCarTolerance)*G4UniformRand(); 
    phi    = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }
  else if ( rand < 2*part )  // Rmin
  {
    radius = pRmin -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 
    pz     = -tubsZ - 0.5*kCarTolerance + (2*tubsZ + kCarTolerance)*G4UniformRand();   
    phi    = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }
  else if ( rand < 3*part )  // phi1
  {
    radius = pRmin - 0.5*kCarTolerance + (pRmax-pRmin+kCarTolerance)*G4UniformRand(); 
    pz     = -tubsZ - 0.5*kCarTolerance + (2*tubsZ + kCarTolerance)*G4UniformRand();   
    phi    = phi1 -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();
  }
  else if ( rand < 4*part )  // phi2
  {
    radius = pRmin - 0.5*kCarTolerance + (pRmax-pRmin+kCarTolerance)*G4UniformRand(); 
    pz     = -tubsZ - 0.5*kCarTolerance + (2*tubsZ + kCarTolerance)*G4UniformRand();   
    phi    = phi2 -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();
  }
  else if ( rand < 5*part )  // -fZ
  {
    radius = pRmin - 0.5*kCarTolerance + (pRmax-pRmin+kCarTolerance)*G4UniformRand(); 

    pz     = -tubsZ - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();   
  
    phi    = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }
  else // fZ
  {
    radius = pRmin - 0.5*kCarTolerance + (pRmax-pRmin+kCarTolerance)*G4UniformRand(); 
    pz     = tubsZ - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();     
    phi    = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }

  px = radius*std::cos(phi);
  py = radius*std::sin(phi);
  
  return G4ThreeVector(px,py,pz);
}

/////////////////////////////////////////////////////////////////////////////
//
// Random vector on cons surface

G4ThreeVector GetVectorOnCons(G4Cons& cons)
{
  G4double phi, pRmin, pRmax, radius, px, py, pz;
  G4double part = 1./6.;
  G4double rand = G4UniformRand();

  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4double kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();

  G4double pRmin1 = cons.GetInnerRadiusMinusZ   ();
  G4double pRmax1 = cons.GetOuterRadiusMinusZ   ();
  G4double pRmin2 = cons.GetInnerRadiusPlusZ   ();
  G4double pRmax2 = cons.GetOuterRadiusPlusZ   ();
  G4double consZ  = cons.GetZHalfLength   ();
  G4double phi1   = cons.GetStartPhiAngle ();
  G4double phi2   = phi1 + cons.GetDeltaPhiAngle ();
  G4double tgMin  = (pRmin2 - pRmin1)/(2.*consZ);
  G4double tgMax  = (pRmax2 - pRmax1)/(2.*consZ);

  if      ( rand < part ) // Rmax
  {
    pz     = -consZ - 0.5*kCarTolerance + (2*consZ + kCarTolerance)*G4UniformRand();
    pRmax  = pRmax1 + tgMax*(pz+consZ); 
    radius = pRmax -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 
    phi    = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }
  else if ( rand < 2*part )  // Rmin
  {
    pz     = -consZ - 0.5*kCarTolerance + (2*consZ + kCarTolerance)*G4UniformRand();   
    pRmin  = pRmin1 + tgMin*(pz+consZ); 
    radius = pRmin -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 
    phi    = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }
  else if ( rand < 3*part )  // phi1
  {
    pz     = -consZ - 0.5*kCarTolerance + (2*consZ + kCarTolerance)*G4UniformRand();   
    pRmax  = pRmax1 + tgMax*(pz+consZ); 
    pRmin  = pRmin1 + tgMin*(pz+consZ); 
    radius = pRmin - 0.5*kCarTolerance + (pRmax-pRmin+kCarTolerance)*G4UniformRand(); 
    phi    = phi1 -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();
  }
  else if ( rand < 4*part )  // phi2
  {
    pz     = -consZ - 0.5*kCarTolerance + (2*consZ + kCarTolerance)*G4UniformRand();   
    pRmax  = pRmax1 + tgMax*(pz+consZ); 
    pRmin  = pRmin1 + tgMin*(pz+consZ); 
    radius = pRmin - 0.5*kCarTolerance + (pRmax-pRmin+kCarTolerance)*G4UniformRand(); 
    phi    = phi2 -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();
  }
  else if ( rand < 5*part )  // -fZ
  {
    pz     = -consZ - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();   
    radius = pRmin1 - 0.5*kCarTolerance + (pRmax1-pRmin1+kCarTolerance)*G4UniformRand();   
    phi    = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }
  else // fZ
  {
    pz     = consZ - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();     
    radius = pRmin2 - 0.5*kCarTolerance + (pRmax2-pRmin2+kCarTolerance)*G4UniformRand(); 
    phi    = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }

  px = radius*std::cos(phi);
  py = radius*std::sin(phi);
  
  return G4ThreeVector(px,py,pz);
}

/////////////////////////////////////////////////////////////////////////////
//
// Random vector on torus surface

G4ThreeVector GetVectorOnTorus(G4Torus& torus)
{
  G4double phi, alpha, radius, px, py, pz;
  G4double part = 1./4.;
  G4double rand = G4UniformRand();

  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4double kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();

  G4double pRmin = torus.GetRmin();
  G4double pRmax = torus.GetRmax();
  G4double pRtor = torus.GetRtor();
  G4double phi1  = torus.GetSPhi();
  G4double phi2  = phi1 + torus.GetDPhi ();

 if      ( rand < part ) // Rmax
  {
    radius = pRmax -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();
    alpha = twopi*G4UniformRand(); 
    phi    = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }
  else if ( rand < 2*part )  // Rmin
  {
    radius = pRmin -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 
    alpha = twopi*G4UniformRand(); 
    phi    = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }
  else if ( rand < 3*part )  // phi1
  {
    radius = pRmin - 0.5*kCarTolerance + (pRmax-pRmin+kCarTolerance)*G4UniformRand(); 
    alpha = twopi*G4UniformRand(); 
    //  pz     = -pRtor - 0.5*kCarTolerance + (2*pRtor + kCarTolerance)*G4UniformRand();   
    phi    = phi1 -0.5*kAngTolerance + (kAngTolerance)*G4UniformRand();
  }
  else   // phi2
  {
    radius = pRmin - 0.5*kCarTolerance + (pRmax-pRmin+kCarTolerance)*G4UniformRand(); 
    alpha = twopi*G4UniformRand(); 
    phi    = phi2 -0.5*kAngTolerance + (kAngTolerance)*G4UniformRand();
  }
  px = (pRtor+radius*std::cos(alpha))*std::cos(phi);
  py = (pRtor+radius*std::cos(alpha))*std::sin(phi);
  pz = radius*std::sin(alpha);

  
  return G4ThreeVector(px,py,pz);
}


enum Esolid {kBox, kOrb, kSphere, kCons, kTubs, kTorus, kPara, kTrapezoid, kTrd};


//////////////////////////////////////////////////////////////////////
//
// Main executable function

int main(int argc, char** argv)
{
  int test_one_solid( Esolid, int, int ); 
  int no_points= 10000;
  int dirs_per_point= 1000; 

  G4cout << "Usage: testSurfaceInOut [ no_surface_points ]  [ no_directions_each ] " 
	 << G4endl << G4endl;
  int points_in=0, dirs_in=0; 
  if( argc >= 2 ){
       points_in = atoi(argv[1]);
  } 
  if( argc >= 3 ){
       dirs_in = atoi(argv[2]);
  }  

  if( points_in > 0 ) { no_points= points_in; } 
  if( dirs_in > 0 ) { dirs_per_point = dirs_in; }  

  G4cout << "Testing each solid with " << no_points << " surface points and " 
	 << dirs_per_point  << " directions each. " << G4endl;

  Esolid useCase; 
  G4cout<< "To test Box." << G4endl;
  test_one_solid( useCase= kBox,  no_points, dirs_per_point ); 

  G4cout<< "To test Tubs." << G4endl;
  test_one_solid( useCase= kTubs,  no_points, dirs_per_point ); 

  G4cout<< "To test Sphere." << G4endl;
  test_one_solid( useCase= kSphere,  no_points, dirs_per_point ); 

  G4cout<< "To test Orb." << G4endl;
  test_one_solid( useCase= kOrb,  no_points, dirs_per_point ); 

  G4cout<< "To test Cons." << G4endl;
  test_one_solid( useCase= kCons,  no_points, dirs_per_point ); 

  G4cout<< "To test Para." << G4endl;
  test_one_solid( useCase= kPara,  no_points, dirs_per_point ); 
 
  G4cout<< "To test Trapezoid." << G4endl;
  test_one_solid( useCase= kTrapezoid,  no_points, dirs_per_point ); 

  G4cout<< "To test Trd." << G4endl;
  test_one_solid( useCase= kTrd,  no_points, dirs_per_point ); 

  return 0; 
}

//

int test_one_solid ( Esolid useCase,  int num_points, int directions_per_point )
{
  G4int i,j;
  G4int iMax=  num_points, jMax= directions_per_point;
  G4int iCheck=iMax/25;
  G4cout << "  Reporting every  " << iCheck <<  " points. " << G4endl; 

  G4double distIn, distOut;

  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  G4double Rtor = 100 ;
  G4double Rmax = Rtor*0.9 ;
  G4double Rmin = Rtor*0.1 ;

  EInside surfaceP;
  G4ThreeVector norm,*pNorm;
  G4bool *pgoodNorm, goodNorm, calcNorm=true;

  pNorm=&norm;
  pgoodNorm=&goodNorm;

  G4cout.precision(20);

  // Boxes for test

  G4Box b1("Test Box #1",20,30,40);
  G4Box b2("Test Box #2",10,10,10);
  G4Box box3("BABAR Box",0.14999999999999999, 
                           24.707000000000001,  
                           22.699999999999999) ;

  // orbs for test

  G4Orb  o1("Solid G4Orb",50);
  G4Orb  o10("s10",0.018*mm);
  // G4Orb* solidO1= new G4Orb("O1", 2.7*cm);


  // spheres for test

  //  G4Sphere s5("Patch (phi/theta seg)",45,50,-pi/4,halfpi,pi/4,halfpi);
  //  G4Sphere s5("Patch (phi/theta seg)",45,50,-pi/4.,pi/4.,pi/4,pi/4.);
  G4Sphere s5("Patch (phi/theta seg)",45,50,-pi/4.,pi/4.,pi/2,pi/4.);

  G4Sphere s6("John example",300,500,0,5.76,0,pi) ; 
  G4Sphere s7("sphere7",1400.,1550.,0.022321428571428572,0.014642857142857141,
	                  1.5631177553663251,0.014642857142857141    );
  G4Sphere s8("sphere",278.746*mm, 280.0*mm, 0.0*degree, 360.0*degree,
		                               0.0*degree, 90.0*degree);

  G4Sphere b216("b216", 1400.0, 1550.0, 
                  0.022321428571428572, 
		  0.014642857142857141,
                  1.578117755366325,
                  0.014642857142867141);

  G4Sphere s9("s9",0*mm,410*mm,0*degree,360*degree,90*degree,90*degree);

  G4Sphere b402("b402", 475*mm, 480*mm, 
               0*degree,360*degree,17.8*degree,144.4*degree);
   

  G4Sphere s10("s10",0*mm,0.018*mm,0*degree,360*degree,0*degree,180*degree);


  G4Sphere s11("s11",5000000.*mm,
                    3700000000.*mm,
                   0*degree,360*degree,0*degree,180*degree);


  G4Sphere sAlex("sAlex",500.*mm,
                    501.*mm,
                   0*degree,360*degree,0*degree,180*degree);

  G4Sphere sLHCB("sLHCB",8600*mm, 8606*mm, 
    -1.699135525184141*degree,
    3.398271050368263*degree,88.52855940538514*degree,2.942881189229715*degree );

  G4Sphere spAroundX("SpAroundX",  10.*mm, 1000.*mm, -1.0*degree, 
                                                     2.0*degree, 
                                        0.*degree, 180.0*degree );


  // G4Tubs t4("Hole Sector #4",45*mm,50*mm,50*mm,halfpi,halfpi);
  G4Tubs t4("Hole Sector #4",45*mm,50*mm,50*mm,-halfpi,halfpi);

  G4Cons c4("Hollow Cut Cone",50,100,50,200,50,-pi/6,pi/3);
  G4Cons c5("Hollow Cut Cone",25,50,75,150,50,0,3*halfpi);
    
  G4Torus torus1("torus1",10.,15.,20.,0.,3*halfpi);
  G4Torus tor2("Hole cutted Torus #2",Rmin,Rmax,Rtor,pi/4,halfpi);
  G4Torus torn2("tn2",Rmin,Rmax,Rtor,halfpi,halfpi);
  G4Torus torn3("tn3",Rmin,Rmax,Rtor,halfpi,3*halfpi);


  // Check of some use cases shown zero-zero problems

  G4ThreeVector pCheck, vCheck;

  // pCheck = G4ThreeVector( 41.418613476008794, -16.893662525384702, 4.9094423552800466 );
  // vCheck = G4ThreeVector( 0.34553222148699703, 0.91172822040596313, 0.22216916084289551 );
  // distIn  = s5.DistanceToIn(pCheck,vCheck);
  // distOut = s5.DistanceToOut(pCheck,vCheck,calcNorm,pgoodNorm,pNorm); 

  // pCheck = G4ThreeVector( 43.169180219772784, -11.564259048580507, 5.2621090605480623 ); 
  // surfaceP = s5.Inside(pCheck);

  pCheck = G4ThreeVector( -51.189087930597751, -17.382942686173514, -45.939946175080983 );
  vCheck = G4ThreeVector( 0.44410267104120671, -0.88563345532911941, -0.13574314117431641 );
  distIn  = c5.DistanceToIn(pCheck,vCheck);
  distOut = c5.DistanceToOut(pCheck,vCheck,calcNorm,pgoodNorm,pNorm); 

  pCheck = G4ThreeVector(-41.407491890396564, -31.155805955816909, -48.18046093035241 );
  vCheck = G4ThreeVector(0.79040557001853884, -0.52472467107317944, -0.31610608100891113 );
  distIn  = c5.DistanceToIn(pCheck,vCheck);
  distOut = c5.DistanceToOut(pCheck,vCheck,calcNorm,pgoodNorm,pNorm); 

  pCheck = G4ThreeVector(-66.68328490196707, -47.460245099793099, -18.151754141035401 );
  vCheck = G4ThreeVector(-0.066679791594195931, 0.88577693677046576, -0.45929622650146484 );
  distIn  = c5.DistanceToIn(pCheck,vCheck);
  distOut = c5.DistanceToOut(pCheck,vCheck,calcNorm,pgoodNorm,pNorm); 


#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif

  
  // Check box tracking function 

  switch (useCase)
  {

    case kBox:
    G4cout<<"Testing G4Box:"<<G4endl<<G4endl;
    for(i=0;i<iMax;i++)
    {
      G4ThreeVector p = GetVectorOnBox(b1);
    
      surfaceP = b1.Inside(p);
      if(surfaceP != kSurface)
      {
        G4cout<<"p is out of surface: "<<G4endl;
        G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl<<G4endl;

      }
      else
      {
        G4ThreeVector origin_norm = b1.SurfaceNormal(p);

        for(j=0;j<jMax;j++)
        {
          G4ThreeVector v = GetRandomUnitVector();

          distIn  = b1.DistanceToIn(p,v);
          distOut = b1.DistanceToOut(p,v,calcNorm,pgoodNorm,pNorm); 

          if( distIn < kCarTolerance && distOut < kCarTolerance )
	  {
	    G4cout<<" distIn < kCarTolerance && distOut < kCarTolerance"<<G4endl;
            G4cout<<"distIn = "<<distIn<<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p: "; // << G4endl;
            G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl;
            G4cout<<" direction v: "; // << G4endl; 
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl<<G4endl;
	  }
          else if(distIn > 100000*kCarTolerance && distOut > 100*kCarTolerance)
	  {
	    G4cout<<" distIn > 100000*kCarTolerance && distOut > 100*kCarTolerance"<<G4endl;
            G4cout<<"distIn = "<<distIn<<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p: ";  // << G4endl;
            G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl;
            G4cout<<" direction v: ";  // << G4endl;
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl;
	    G4cout<<" dot prod v.norm(orig)= " << v.dot( origin_norm ) << G4endl;
	    G4cout <<G4endl;
	  }     
        }
      }
    }
    break;
    
    case kOrb:
      G4cout<<"Testing G4Orb:"<<G4endl<<G4endl;
    for(i=0;i<iMax;i++)
    {
      G4ThreeVector p = GetVectorOnOrb(o1);
    
      surfaceP = o1.Inside(p);
      if(surfaceP != kSurface)
      {
        G4cout<<"p is out of surface: "<<G4endl;
        G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl<<G4endl;

      }
      else
      {
        for(j=0;j<jMax;j++)
        {
          G4ThreeVector v = GetRandomUnitVector();

          distIn  = o1.DistanceToIn(p,v);
          distOut = o1.DistanceToOut(p,v,calcNorm,pgoodNorm,pNorm); 

          if( distIn < kCarTolerance && distOut < kCarTolerance )
	  {
	    G4cout<<" distIn < kCarTolerance && distOut < kCarTolerance"<<G4endl;
            G4cout<<"distIn = "<<distIn<<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p: "<<G4endl;
            G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl;
            G4cout<<" direction v: "<<G4endl;
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl<<G4endl;
	  }
          else if(distIn > 100000*kCarTolerance && distOut > 100*kCarTolerance)
	  {
	    G4cout<<" distIn > 100000*kCarTolerance && distOut > 100*kCarTolerance"<<G4endl;
            G4cout<<"distIn = "<<distIn<<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p: "<<G4endl;
            G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl;
            G4cout<<" direction v: "<<G4endl;
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl; 
	    G4cout<<" p.v/|p|= " << p.dot(v)/p.mag() <<G4endl;
	  }     
        }
      }
    }
    break;

    case kSphere:
      G4cout<<"Testing one cutted G4Sphere:" << s5 <<G4endl<<G4endl;
    for(i=0;i<iMax;i++)
    {
      if(i%iCheck == 0) G4cout<<"i = "<<i<<G4endl;
      G4ThreeVector p = GetVectorOnSphere(s5);      
      surfaceP = s5.Inside(p);
      if(surfaceP != kSurface)
      {
        G4cout<<"p is out of surface: "<<G4endl;
        G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl<<G4endl;

      }
      else
      {
        for(j=0;j<jMax;j++)
        {
          G4ThreeVector v = GetRandomUnitVector();

          distIn  = s5.DistanceToIn(p,v);
          distOut = s5.DistanceToOut(p,v,calcNorm,pgoodNorm,pNorm); 

	  //  if( distIn < kCarTolerance && distOut < kCarTolerance )
          if( distIn == 0. && distOut == 0. )
	  {
	    G4cout<<" distIn < kCarTolerance && distOut < kCarTolerance"<<G4endl;
            G4cout<<"distIn = "<<distIn<<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p: "<<G4endl;
            G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl;
            G4cout<<" direction v: "<<G4endl;
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl<<G4endl;
	  }
          else if(distIn > 100000*kCarTolerance && distOut > 100*kCarTolerance)
	  {
	    /*
	    G4cout<<" distIn > 100000*kCarTolerance && distOut > 100*kCarTolerance"<<G4endl;
            G4cout<<"distIn = "<<distIn<<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p: "<<G4endl;
            G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl;
            G4cout<<" direction v: "<<G4endl;
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl<<G4endl;
	    */
	  }     
        }
      }
    }
    break;

    case kTubs:
      G4cout<<"Testing all cutted G4Tubs:"<<G4endl<<G4endl;
    for(i=0;i<iMax;i++)
    {
      if(i%iCheck == 0) G4cout<<"i = "<<i<<G4endl;
      G4ThreeVector p = GetVectorOnTubs(t4);      
      surfaceP = t4.Inside(p);
      if(surfaceP != kSurface)
      {
        G4cout<<"p is out of surface: "<<G4endl;
        G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl<<G4endl;

      }
      else
      {
        for(j=0;j<jMax;j++)
        {
          G4ThreeVector v = GetRandomUnitVector();

          distIn  = t4.DistanceToIn(p,v);
          distOut = t4.DistanceToOut(p,v,calcNorm,pgoodNorm,pNorm); 

	  //  if( distIn < kCarTolerance && distOut < kCarTolerance )
          if( distIn == 0. && distOut == 0. )
	  {
	    G4cout<<" distIn < kCarTolerance && distOut < kCarTolerance"<<G4endl;
            G4cout<<"distIn = "<<distIn<<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p: "<<G4endl;
            G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl;
            G4cout<<" direction v: "<<G4endl;
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl<<G4endl;
	  }
          else if(distIn > 100000*kCarTolerance && distOut > 100*kCarTolerance)
	  {
	    /*
	    G4cout<<" distIn > 100000*kCarTolerance && distOut > 100*kCarTolerance"<<G4endl;
            G4cout<<"distIn = "<<distIn<<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p: "<<G4endl;
            G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl;
            G4cout<<" direction v: "<<G4endl;
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl<<G4endl;
	    */
	  }     
        }
      }
    }
    break;

    case kCons:
      G4cout<<"Testing all cutted G4Cons:"<<G4endl<<G4endl;
    for(i=0;i<iMax;i++)
    {
      if(i%iCheck == 0) G4cout<<"i = "<<i<<G4endl;
      G4ThreeVector p = GetVectorOnCons(c5);      
      surfaceP = c5.Inside(p);
      if(surfaceP != kSurface)
      {
        G4cout<<"p is out of surface: "<<G4endl;
        G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl<<G4endl;

      }
      else
      {
        for(j=0;j<jMax;j++)
        {
          G4ThreeVector v = GetRandomUnitVector();

          distIn  = c5.DistanceToIn(p,v);
          distOut = c5.DistanceToOut(p,v,calcNorm,pgoodNorm,pNorm); 

	  //  if( distIn < kCarTolerance && distOut < kCarTolerance )
          if( distIn == 0. && distOut == 0. )
	  {
	    G4cout<<" distIn < kCarTolerance && distOut < kCarTolerance"<<G4endl;
            G4cout<<"distIn = "<<distIn<<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p: "<<G4endl;
            G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl;
            G4cout<<" direction v: "<<G4endl;
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl<<G4endl;
	  }
          else if(distIn > 100000*kCarTolerance && distOut > 100*kCarTolerance)
	  {
	    /*
	    G4cout<<" distIn > 100000*kCarTolerance && distOut > 100*kCarTolerance"<<G4endl;
            G4cout<<"distIn = "<<distIn<<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p: "<<G4endl;
            G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl;
            G4cout<<" direction v: "<<G4endl;
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl<<G4endl;
	    */
	  }     
        }
      }
    }
    break;
    case kTorus:
      G4cout<<"Testing all cutted G4Torus:"<<G4endl<<G4endl;
    for(i=0;i<iMax;i++)
    {
      if(i%iCheck == 0) G4cout<<"i = "<<i<<G4endl;
      G4ThreeVector p = GetVectorOnTorus(torus1);      
      surfaceP = torus1.Inside(p);
      if(surfaceP != kSurface)
      {
        G4cout<<"p is out of surface: "<<G4endl;
        G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl<<G4endl;

      }
      else
      {
        for(j=0;j<jMax;j++)
        {
          G4ThreeVector v = GetRandomUnitVector();

          distIn  = torus1.DistanceToIn(p,v);
          distOut = torus1.DistanceToOut(p,v,calcNorm,pgoodNorm,pNorm); 

	  //  if( distIn < kCarTolerance && distOut < kCarTolerance )
          if( distIn == 0. && distOut == 0. )
	  {
	    G4cout<<" distIn < kCarTolerance && distOut < kCarTolerance"<<G4endl;
            G4cout<<"distIn = "<<distIn<<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p: "<<G4endl;
            G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl;
            G4cout<<" direction v: "<<G4endl;
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl<<G4endl;
	  }
          else if(distIn > 100000*kCarTolerance && distOut > 100*kCarTolerance)
	  {
	    /*
	    G4cout<<" distIn > 100000*kCarTolerance && distOut > 100*kCarTolerance"<<G4endl;
            G4cout<<"distIn = "<<distIn<<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p: "<<G4endl;
            G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl;
            G4cout<<" direction v: "<<G4endl;
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl<<G4endl;
	    */
	  }     
        }
      }
    }
    break;

    default:
      G4cout<<"A test is not implemented for this case: " << useCase <<G4endl<<G4endl;

    break;   
  }
  return 0;
}

