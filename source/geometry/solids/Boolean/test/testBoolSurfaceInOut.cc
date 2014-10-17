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
// 14.07.04 V.Grichine creation for box, orb and sphere
// 15.02.05 V.Grichine changes for boolean solids

#include "G4ios.hh"
#include <assert.h>
#include <cmath>
#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "geomdefs.hh"
#include "Randomize.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"
#include "G4GeometryTolerance.hh"

#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
// #include "G4Hype.hh"
#include "G4Para.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"

#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"


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

  a = box.GetXHalfLength();
  b = box.GetYHalfLength();
  c = box.GetZHalfLength();
  
  rand = G4UniformRand();
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

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

  G4double pRmin  = sphere.GetInnerRadius();
  G4double pRmax  = sphere.GetOuterRadius();
  G4double phi1   = sphere.GetStartPhiAngle();
  G4double phi2   = phi1 + sphere.GetDeltaPhiAngle();
  G4double theta1 = sphere.GetStartThetaAngle();
  G4double theta2 = theta1 + sphere.GetDeltaThetaAngle();
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4double kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();

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

  G4double pRmin = tubs.GetInnerRadius   ();
  G4double pRmax = tubs.GetOuterRadius   ();
  G4double tubsZ = tubs.GetZHalfLength   ();
  G4double phi1  = tubs.GetStartPhiAngle ();
  G4double phi2  = phi1 + tubs.GetDeltaPhiAngle ();
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4double kAngTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

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
// Random vector out of tubs surface

G4ThreeVector GetVectorOutOfTubs(G4Tubs& tubs)
{
  G4double phi, radius, px, py, pz;
  G4double part = 1./5.;
  G4double rand = G4UniformRand();

  G4double pRmin = tubs.GetInnerRadius   ();
  G4double pRmax = tubs.GetOuterRadius   ();
  G4double tubsZ = tubs.GetZHalfLength   ();
  G4double phi1  = tubs.GetStartPhiAngle ();
  G4double deltaPhi  = tubs.GetDeltaPhiAngle ();
  G4double phi2  = phi1 + deltaPhi;
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  if      ( rand < part ) // Rmax
  {
    radius = pRmax +0.5*kCarTolerance + pRmax*G4UniformRand(); 
    pz     = -2*tubsZ + (4*tubsZ)*G4UniformRand(); 
    phi    = (2*pi)*G4UniformRand();
  }
  else if ( rand < 2*part )  // Rmin
  {
    radius = (pRmin -0.5*kCarTolerance)*G4UniformRand(); 
    pz     = -2*tubsZ + (4*tubsZ)*G4UniformRand(); 
    phi    = (2*pi)*G4UniformRand();
  }
  else if ( rand < 3*part )  // out of deltaPhi
  {
    radius = (2*pRmax)*G4UniformRand(); 
    pz     = -2*tubsZ - 0.5*kCarTolerance + (4*tubsZ + kCarTolerance)*G4UniformRand();   
    phi    = phi2 + 0.5*kCarTolerance + (2*pi - deltaPhi - kCarTolerance)*G4UniformRand();
  }
  else if ( rand < 4*part )  // -fZ
  {
    radius = (2*pRmax)*G4UniformRand(); 
    pz     = -tubsZ - 0.5*kCarTolerance - (tubsZ)*G4UniformRand();   
    phi    = (2*pi)*G4UniformRand();  
  }
  else // fZ
  {
    radius = (2*pRmax)*G4UniformRand(); 
    pz     = tubsZ + 0.5*kCarTolerance + (tubsZ)*G4UniformRand();   
    phi    = (2*pi)*G4UniformRand();  
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

  G4double pRmin1 = cons.GetInnerRadiusMinusZ   ();
  G4double pRmax1 = cons.GetOuterRadiusMinusZ   ();
  G4double pRmin2 = cons.GetInnerRadiusPlusZ   ();
  G4double pRmax2 = cons.GetOuterRadiusPlusZ   ();
  G4double consZ  = cons.GetZHalfLength   ();
  G4double phi1   = cons.GetStartPhiAngle ();
  G4double phi2   = phi1 + cons.GetDeltaPhiAngle ();
  G4double tgMin  = (pRmin2 - pRmin1)/(2.*consZ);
  G4double tgMax  = (pRmax2 - pRmax1)/(2.*consZ);
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4double kAngTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

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
  G4double phi, radius, px, py, pz;
  G4double part = 1./4.;
  G4double rand = G4UniformRand();

  G4double pRmin = torus.GetRmin();
  G4double pRmax = torus.GetRmax();
  G4double pRtor = torus.GetRtor();
  G4double phi1  = torus.GetSPhi();
  G4double phi2  = phi1 + torus.GetDPhi ();
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4double kAngTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  if      ( rand < part ) // Rmax
  {
    radius = pRmax -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 
    pz     = -pRtor - 0.5*kCarTolerance + (2*pRtor + kCarTolerance)*G4UniformRand(); 
    phi    = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }
  else if ( rand < 2*part )  // Rmin
  {
    radius = pRmin -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 
    pz     = -pRtor - 0.5*kCarTolerance + (2*pRtor + kCarTolerance)*G4UniformRand();   
    phi    = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }
  else if ( rand < 3*part )  // phi1
  {
    radius = pRmin - 0.5*kCarTolerance + (pRmax-pRmin+kCarTolerance)*G4UniformRand(); 
    pz     = -pRtor - 0.5*kCarTolerance + (2*pRtor + kCarTolerance)*G4UniformRand();   
    phi    = phi1 -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();
  }
  else if ( rand < 4*part )  // phi2
  {
    radius = pRmin - 0.5*kCarTolerance + (pRmax-pRmin+kCarTolerance)*G4UniformRand(); 
    pz     = -pRtor - 0.5*kCarTolerance + (2*pRtor + kCarTolerance)*G4UniformRand();   
    phi    = phi2 -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();
  }
  else if ( rand < 5*part )  // -fZ
  {
    radius = pRmin - 0.5*kCarTolerance + (pRmax-pRmin+kCarTolerance)*G4UniformRand(); 

    pz     = -pRtor - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();   
  
    phi    = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }
  else // fZ
  {
    radius = pRmin - 0.5*kCarTolerance + (pRmax-pRmin+kCarTolerance)*G4UniformRand(); 
    pz     = pRtor - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();     
    phi    = phi1 - 0.5*kAngTolerance + (phi2 - phi1 + kAngTolerance)*G4UniformRand();
  }

  px = radius*std::cos(phi);
  py = radius*std::sin(phi);
  
  return G4ThreeVector(px,py,pz);
}


//////////////////////////////////////////////////////////////////////
//
// Main executable function

int main(void)
{
  G4int i, j, iMax=1000, jMax=1000;
  G4int iCheck = iMax/10;
  G4double distIn, distIn1, distIn2, distOut;
  EInside surfaceP, surfaceP1, surfaceP2;
  G4ThreeVector norm, *pNorm;
  G4bool *pgoodNorm, goodNorm, calcNorm=true;
  G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  enum Esolid {kInter, kBox, kOrb, kSphere, kCons, kTubs, kTorus, kPara, kTrapezoid, kTrd};

  Esolid useCase = kInter;

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
 


    G4Tubs*  detTub1 = new G4Tubs("Tubs4", 0.*cm, 5.*cm, 5.*cm, 0., 359*deg);
    G4Tubs*  detTub12 = new G4Tubs("Tubs12", 1.*cm, 5.*cm, 5.*cm, 0., 359*deg);
    G4Tubs*  detTub2 = new G4Tubs("Tubs5", 1.*cm, 5.1*cm, 5.1*cm, 0., 359*deg);

    G4IntersectionSolid*  detInt1 = new G4IntersectionSolid("inter1", 
                                        detTub1, detTub2, 0, G4ThreeVector() );
    G4IntersectionSolid*  detInt2 = new G4IntersectionSolid("inter2", 
                                        detTub2, detTub1, 0, G4ThreeVector() );

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

  pCheck = G4ThreeVector( -17.140591059524617, -23.320101452294466, -49.999999999668375 ); 
  vCheck = G4ThreeVector ( -0.69080640316788089, -0.58982688856527554, 0.41819941997528076 ); 
  distIn  = detTub12->DistanceToIn(pCheck,vCheck);
  distIn1  = detTub1->DistanceToIn(pCheck,vCheck);
  distIn2  = detTub2->DistanceToIn(pCheck,vCheck);
  G4cout<<"dTub12 = "<<distIn<<";    dTub1 = "<<distIn1<<";    dTub2 = "<<distIn2<<G4endl;
  surfaceP   = detTub12->Inside(pCheck);
  surfaceP1  = detTub1->Inside(pCheck);
  surfaceP2  = detTub2->Inside(pCheck);
  G4cout<<"insideTub12 = "<<OutputInside(surfaceP)<<";     insideTub1 = "<<OutputInside(surfaceP1)
        <<";     insideTub2 = "<<OutputInside(surfaceP2)<<G4endl; 
  distIn1  = detInt1->DistanceToIn(pCheck,vCheck);
  distIn2  = detInt2->DistanceToIn(pCheck,vCheck);
  G4cout<<"Int1 = "<<distIn1<<";        Int2 = "<<distIn2<<G4endl;
  surfaceP1  = detInt1->Inside(pCheck);
  surfaceP2  = detInt2->Inside(pCheck);
  G4cout<<"insideInt1 = "<<OutputInside(surfaceP1)
        <<";       insideInt2 = "<<OutputInside(surfaceP2)<<G4endl; 


#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
  
  // Check box tracking function 

  switch (useCase)
  {
    case kInter:
    G4cout<<"Testing of all cutted G4Tubs intersection:"<<G4endl<<G4endl;
    for( i = 0; i < iMax; i++ )
    {
      if(i%iCheck == 0) G4cout<<"i = "<<i<<G4endl;
          
      // G4ThreeVector p1 = GetVectorOnTubs(*detTub1);
      G4ThreeVector p1 = GetVectorOutOfTubs(*detTub12);
      G4ThreeVector p2 = GetVectorOnTubs(*detTub2);

      surfaceP = detInt1->Inside(p1);

      if(
         //  surfaceP != kSurface
           surfaceP != kOutside
                                  )
      {
        // G4cout<<"p is out of surface: "<<G4endl;
        // G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl<<G4endl;

      }
      else
      {
        for( j = 0; j < jMax; j++ )
        {
          G4ThreeVector v = GetRandomUnitVector();

          distIn1  = detInt1->DistanceToIn(p1,v);
          distIn2  = detInt2->DistanceToIn(p1,v);
          // distOut  = detInt1->DistanceToOut(p1,v,calcNorm,pgoodNorm,pNorm);
          // distOut = t4.DistanceToOut(p,v,calcNorm,pgoodNorm,pNorm); 

	  //  if( distIn < kCarTolerance && distOut < kCarTolerance )
          if( 
	     // distIn1 !=  distIn2
	     std::abs( distIn1 -  distIn2 ) > 100*kCarTolerance
              //   distIn1 ==  distOut 
                                          )
	  {
	    G4cout<<" distIn1  != distIn2: "<<G4endl;
	    //   G4cout<<" distIn1  == distOut: "<<G4endl;
            G4cout<<"distIn1 = "<<distIn1
	          <<";  distIn2 = "<<distIn2<<G4endl;
	      //  <<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p1: "<<G4endl;
            G4cout<<"( "<<p1.x()<<", "<<p1.y()<<", "<<p1.z()<<" ); "<<G4endl;
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
        for(j=0;j<jMax;j++)
        {
          G4ThreeVector v = GetRandomUnitVector();

          distIn  = b1.DistanceToIn(p,v);
          distOut = b1.DistanceToOut(p,v,calcNorm,pgoodNorm,pNorm); 

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
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl<<G4endl;
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
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl<<G4endl;
	  }     
        }
      }
    }
    break;

    case kSphere:
      G4cout<<"Testing all cutted G4Sphere:"<<G4endl<<G4endl;
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
    break;   
  }
  return 0;
}

