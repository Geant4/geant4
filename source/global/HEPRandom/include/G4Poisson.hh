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
// $Id: G4Poisson.hh,v 1.6 2001-07-11 10:00:43 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
// ------------------------------------------------------------
// Class description:
//
// G4Poisson is the C++ implementation of the CERNLIB GPOISS algorithm
// for the generation of Poisson distributed random numbers. It has been
// adapted to invoke HepRandom from CLHEP for the primary engine generators.
// GPOISS is recognized to be a faster algorithm, providing however a less
// accurate output, than the algorithm adopted in CLHEP.

// ------------------------------------------------------------
#ifndef G4POISSON_HH
#define G4POISSON_HH

#include "globals.hh"
#include "Randomize.hh"

inline G4long G4Poisson(G4double mean)
{
  G4long number = 0;
  const G4int border = 16;
  G4double limit = 2e9;

  if(mean <= border) {
    G4double position = RandFlat::shoot();
    G4double poissonValue = exp(-mean);
    G4double poissonSum = poissonValue;

    while(poissonSum <= position) {
      number++ ;
      poissonValue *= mean/number;
      poissonSum += poissonValue;
    }
    return number;
  } // the case of mean <= 16

  G4double value, t, y;
  t = sqrt(-2*log(RandFlat::shoot()));
  y = twopi*RandFlat::shoot();
  t *= cos(y);
  value = mean + t*sqrt(mean) + 0.5;
  if(value <= 0) {return 0;}
  if(value >= limit) { return G4long(limit);}
  return G4long(value);
}

#endif  /* G4POISSON_HH */
