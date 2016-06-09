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
// 07.11.08 - P.Gumplinger, based on implementation in G4OpBoundaryProcess
//
// ---------------------------------------------------------------------------

#ifndef G4RANDOMTOOLS_HH
#define G4RANDOMTOOLS_HH

#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "Randomize.hh"
#include "G4RandomDirection.hh"

// ---------------------------------------------------------------------------
// Returns a random lambertian unit vector
//
inline G4ThreeVector G4LambertianRand(const G4ThreeVector& normal)
{
  G4ThreeVector vect;
  G4double ndotv;

  do
  {
    vect = G4RandomDirection();
    ndotv = normal * vect;

    if (ndotv < 0.0)
    {
      vect = -vect;
      ndotv = -ndotv;
    }

  } while (!(G4UniformRand() < ndotv));

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

#endif  /* G4RANDOMTOOLS_HH */
