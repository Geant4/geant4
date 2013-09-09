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
// $Id:$

// ------------------------------------------------------------

#include <CLHEP/Units/PhysicalConstants.h>

#include "G4Poisson.hh"
#include "G4Exp.hh"
#include "Randomize.hh"

G4long G4Poisson(G4double mean)
{
  G4long number = 0;
  const G4int border = 16;
  G4double limit = 2e9;

  if(mean <= border)
  {
    G4double position = G4UniformRand();
    G4double poissonValue = G4Exp(-mean);
    G4double poissonSum = poissonValue;

    while(poissonSum <= position)
    {
      number++ ;
      poissonValue *= mean/number;
      poissonSum += poissonValue;
    }
    return number;
  } // the case of mean <= 16

  G4double value, t, y;
  t = std::sqrt(-2*std::log(G4UniformRand()));
  y = CLHEP::twopi*G4UniformRand();
  t *= std::cos(y);
  value = mean + t*std::sqrt(mean) + 0.5;
  if(value <= 0) {return 0;}
  if(value >= limit) { return G4long(limit);}
  return G4long(value);
}
