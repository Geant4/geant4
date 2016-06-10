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
// $Id: G4Poisson.hh 79932 2014-03-26 10:43:05Z gcosmo $
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

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4Types.hh"
#include "G4Exp.hh"
#include "Randomize.hh"

inline G4long G4Poisson(G4double mean)
{
  G4long number = 0;
  const G4int border = 16;
  const G4double limit = 2e9;

  if(mean <= border)
  {
    G4double position = G4UniformRand();
    G4double poissonValue = G4Exp(-mean);
    G4double poissonSum = poissonValue;

    while(poissonSum <= position)
    {
      ++number;
      poissonValue *= mean/number;
      poissonSum += poissonValue;
    }
    return number;
  } // the case of mean <= 16

  G4double t = std::sqrt(-2.*std::log(G4UniformRand()))*
               std::cos(2.*CLHEP::pi*G4UniformRand());
  G4double value = mean + t*std::sqrt(mean) + 0.5;
  if(value < 0.) {return 0;}
  return (value >= limit) ? G4long(limit) : G4long(value);
}

#endif  /* G4POISSON_HH */
