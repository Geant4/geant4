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
// $Id$
//
// Abstract base class for multibody uniform phase space generators.
// Subclasses implement a specific algorithm, such as Kopylov, GENBOD,
// or Makoto's NBody.  Subclasses are used by G4HadDecayGenerator.
//
// Author:	Michael Kelsey (SLAC) <kelsey@slac.stanford.edu>

#include "G4VHadPhaseSpaceAlgorithm.hh"
#include "G4HadronicException.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <vector>




// Two body decay with uniform angular distribution

void G4VHadPhaseSpaceAlgorithm::
GenerateTwoBody(G4double initialMass,
		const std::vector<G4double>& masses,
		std::vector<G4LorentzVector>& finalState) {
  if (GetVerboseLevel()>1) 
    G4cout << " >>> G4HadDecayGenerator::FillTwoBody" << G4endl;

  // Initialization and sanity check
  finalState.clear();
  if (masses.size() != 2U) return;	// Should not have been called

  // Momentum of final state (energy balance has already been checked)
  G4double p = TwoBodyMomentum(initialMass,masses[0],masses[1]);
  if (GetVerboseLevel()>2) G4cout << " finalState momentum = " << p << G4endl;

  finalState.resize(2);				// Allows filling by index
  finalState[0].setVectM(UniformVector(p), masses[0]);
  finalState[1].setVectM(-finalState[0].vect(), masses[1]);
}


// Samples a random vector with given magnitude

G4ThreeVector G4VHadPhaseSpaceAlgorithm::UniformVector(G4double mag) const {
  // FIXME:  Should this be made a static thread-local buffer?
  G4ThreeVector v;
  v.setRThetaPhi(mag, UniformTheta(), UniformPhi());
  return v;
}
