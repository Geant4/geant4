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
// $Id: G4RandomTools.hh 105368 2017-07-24 09:44:20Z gcosmo $
//
// 
// ---------------------------------------------------------------------------
//      GEANT 4 class header file 
// ---------------------------------------------------------------------------
// Class description:
//
// Utility functions 

// History:
//
// 24.08.17 - E.Tcherniaev, added G4RandomRadiusInRing, G4RandomPointInEllipse
//                          G4RandomPointOnEllipse, G4RandomPointOnEllipsoid
// 07.11.08 - P.Gumplinger, based on implementation in G4OpBoundaryProcess
//
// ---------------------------------------------------------------------------

#ifndef G4RANDOMTOOLS_HH
#define G4RANDOMTOOLS_HH

#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "Randomize.hh"
#include "G4TwoVector.hh"
#include "G4ThreeVector.hh"
#include "G4RandomDirection.hh"

// ---------------------------------------------------------------------------
// Returns a random lambertian unit vector (rejection sampling)
//
inline G4ThreeVector G4LambertianRand(const G4ThreeVector& normal)
{
  G4ThreeVector vect;
  G4double ndotv;
  G4int count=0;
  const G4int max_trials = 1024;

  do
  {
    ++count;
    vect = G4RandomDirection();
    ndotv = normal * vect;

    if (ndotv < 0.0)
    {
      vect = -vect;
      ndotv = -ndotv;
    }

  } while (!(G4UniformRand() < ndotv) && (count < max_trials));

  return vect;
}

// ---------------------------------------------------------------------------
// Chooses a random vector within a plane given by the unit normal
//
inline G4ThreeVector G4PlaneVectorRand(const G4ThreeVector& normal)
{
  G4ThreeVector vec1 = normal.orthogonal();
  G4ThreeVector vec2 = vec1.cross(normal);

  G4double phi = CLHEP::twopi*G4UniformRand();
  G4double cosphi = std::cos(phi);
  G4double sinphi = std::sin(phi);

  return cosphi * vec1 + sinphi * vec2;
}

// ---------------------------------------------------------------------------
// Returns a random radius in annular ring
//
inline G4double G4RandomRadiusInRing(G4double rmin, G4double rmax)
{
  if (rmin == rmax)
  {
    return rmin;
  }
  G4double k = G4UniformRand();
  return (rmin <= 0) ? rmax*std::sqrt(k)
                     : std::sqrt(k*rmax*rmax + (1.-k)*rmin*rmin);
}

// ---------------------------------------------------------------------------
// Returns a random point in ellipse (x/a)^2 + (y/b)^2 = 1
// (rejection sampling)
//
inline G4TwoVector G4RandomPointInEllipse(G4double a, G4double b)
{
  G4double aa = (a*a == 0) ? 0 : 1/(a*a);
  G4double bb = (b*b == 0) ? 0 : 1/(b*b);
  for (G4int i=0; i<1000; ++i)
  {
    G4double x = a*(2*G4UniformRand() - 1);
    G4double y = b*(2*G4UniformRand() - 1);
    if (x*x*aa + y*y*bb <= 1) return G4TwoVector(x,y);
  }
  return G4TwoVector(0,0);
}

// ---------------------------------------------------------------------------
// Returns a random point on ellipse (x/a)^2 + (y/b)^2 = 1
// (rejection sampling)
//
inline G4TwoVector G4RandomPointOnEllipse(G4double a, G4double b)
{
  G4double A = std::abs(a);
  G4double B = std::abs(b);
  G4double mu_max = std::max(A,B);

  G4double x,y;
  for (G4int i=0; i<1000; ++i)
  {
    G4double phi = CLHEP::twopi*G4UniformRand();
    x = std::cos(phi);
    y = std::sin(phi);
    G4double mu = std::sqrt((B*x)*(B*x) + (A*y)*(A*y));
    if (mu_max*G4UniformRand() <= mu) break;
  }
  return G4TwoVector(A*x,B*y);
}

// ---------------------------------------------------------------------------
// Returns a random point on ellipsoid (x/a)^2 + (y/b)^2 + (z/c)^2 = 1
// (rejection sampling)
//
inline
G4ThreeVector G4RandomPointOnEllipsoid(G4double a, G4double b, G4double c)
{
  G4double A = std::abs(a);
  G4double B = std::abs(b);
  G4double C = std::abs(c);
  G4double mu_max = std::max(std::max(A*B,A*C),B*C);

  G4ThreeVector p;
  for (G4int i=0; i<1000; ++i)
  {
    p = G4RandomDirection();
    G4double xbc = p.x()*B*C;
    G4double yac = p.y()*A*C;
    G4double zab = p.z()*A*B;
    G4double mu = std::sqrt(xbc*xbc + yac*yac + zab*zab);
    if (mu_max*G4UniformRand() <= mu) break;
  }
  return G4ThreeVector(A*p.x(),B*p.y(),C*p.z());
}

#endif  /* G4RANDOMTOOLS_HH */
