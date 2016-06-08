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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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

G4DummyProbability::G4DummyProbability(const G4DummyProbability& right)
{

  G4Exception("G4DummyProbability::copy_constructor meant to not be accessible");

}

const G4DummyProbability& G4DummyProbability::
operator=(const G4DummyProbability& right) 
{

  G4Exception("G4DummyProbability::operator= meant to not be accessible");
  return *this;
}

G4bool G4DummyProbability::operator==(const G4DummyProbability& right) const
{

  return false;

}

G4bool G4DummyProbability::operator!=(const G4DummyProbability& right) const
{

  return true;

}

// Calculate the emission probability
//

G4double G4DummyProbability::EmissionProbDensity(const G4Fragment& frag, 
                                                 const G4double exciteE)
{

  G4double theProb = 0.0;

  return theProb;

}

G4double G4DummyProbability::EmissionProbability(const G4Fragment& frag, 
                                                 const G4double exciteE)
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


