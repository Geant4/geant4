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
// $Id: G4RandomDirection.hh,v 1.2 2004/11/12 17:38:34 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
// ------------------------------------------------------------
// Class description:
//
// Funtion returning a unit 3-vector homogeneously randomised over 4pi
// solid angle. It can be used in any particle scattering methods
// instead of:
//   z=R1, x=SQRT(1-R1*R1)*SIN(2*pi*R2), y=SQRT(1-R1*R1)*COS(2*pi*R2)
// providing more performant results.
// ------------------------------------------------------------
#ifndef G4RANDOMDIR_HH
#define G4RANDOMDIR_HH

#include "globals.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"

inline G4ThreeVector G4RandomDirection()
{
  // Randomization in one of 8 Quadrants (x>0, y>0, z>0)
  //
  G4double x=G4UniformRand(), y=G4UniformRand(), z=G4UniformRand();
  G4double r2= x*x+y*y+z*z;
  while(r2>1.||r2<.000001)
  {
    x = G4UniformRand(); y = G4UniformRand(); z = G4UniformRand();
    r2=x*x+y*y+z*z;
  }
  G4double r=std::sqrt(r2), quad=G4UniformRand();

  if(quad>0.5)
  {
    if(quad>0.75)
    {
      if(quad>0.875)    return G4ThreeVector(-x/r,-y/r,-z/r);
      else              return G4ThreeVector(-x/r,-y/r, z/r);
    }
    else
    {
      if(quad>0.625)    return G4ThreeVector(-x/r, y/r,-z/r);
      else              return G4ThreeVector(-x/r, y/r, z/r);
    }
  }
  else
  {
    if(quad>0.25)
    {
      if(quad>0.375)    return G4ThreeVector( x/r,-y/r,-z/r);
      else              return G4ThreeVector( x/r,-y/r, z/r);
    }
    else if(quad>0.125) return G4ThreeVector( x/r, y/r,-z/r);
  }
  return                       G4ThreeVector( x/r, y/r, z/r);
}

#endif  /* G4RANDOMDIR_HH */
