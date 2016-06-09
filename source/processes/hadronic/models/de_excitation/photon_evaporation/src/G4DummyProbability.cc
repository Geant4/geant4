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
//  Class G4DummyProbability.cc
//

#include "G4DummyProbability.hh"
#include "G4ConstantLevelDensityParameter.hh"
#include "Randomize.hh"

// Constructors and operators
//

G4DummyProbability::G4DummyProbability(const G4DummyProbability& ) : G4VEmissionProbability()
{

  throw G4HadronicException(__FILE__, __LINE__, "G4DummyProbability::copy_constructor meant to not be accessible");

}

const G4DummyProbability& G4DummyProbability::
operator=(const G4DummyProbability& ) 
{

  throw G4HadronicException(__FILE__, __LINE__, "G4DummyProbability::operator= meant to not be accessible");
  return *this;
}

G4bool G4DummyProbability::operator==(const G4DummyProbability& ) const
{

  return false;

}

G4bool G4DummyProbability::operator!=(const G4DummyProbability& ) const
{

  return true;

}

// Calculate the emission probability
//

G4double G4DummyProbability::EmissionProbDensity(const G4Fragment& , 
                                                 const G4double )
{

  G4double theProb = 0.0;

  return theProb;

}

G4double G4DummyProbability::EmissionProbability(const G4Fragment& , 
                                                 const G4double )
{

  // From nuclear fragment properties and the excitation energy, calculate
  // the probability for photon evaporation down to last ground level.
  // fragment = nuclear fragment BEFORE de-excitation

  G4double theProb = 0.0;

  // Fall-back is a uniform random number

  G4double uniformNum = G4UniformRand();
  theProb = uniformNum;

  return theProb;

}

G4DummyProbability::~G4DummyProbability() {}


