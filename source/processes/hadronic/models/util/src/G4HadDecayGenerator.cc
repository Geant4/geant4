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
// Multibody "phase space" generator, which provides multiple algorithms
// for sampling.  Momentum vectors are generated in the center-of-mass
// frame of the decay, and returned in a user-supplied buffer.  A sampling
// algorithm is specified via constructor argument.
//
// Author:	Michael Kelsey (SLAC) <kelsey@slac.stanford.edu>

#include "G4HadDecayGenerator.hh"
#include "G4VHadDecayAlgorithm.hh"
#include "G4HadPhaseSpaceKopylov.hh"
#include "G4HadPhaseSpaceGenbod.hh"
#include "G4HadPhaseSpaceNBodyAsai.hh"
#include "G4HadronicException.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <vector>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <iostream>


// Constructors and destructor

G4HadDecayGenerator::G4HadDecayGenerator(Algorithm alg, G4int verbose)
  : verboseLevel(verbose), theAlgorithm(0) {
  switch (alg) {
  case Kopylov: theAlgorithm = new G4HadPhaseSpaceKopylov(verboseLevel); break;
  case GENBOD: theAlgorithm = new G4HadPhaseSpaceGenbod(verboseLevel); break;
  case NBody: theAlgorithm = new G4HadPhaseSpaceNBodyAsai(verboseLevel); break;
  case NONE: theAlgorithm = 0; break;	// User may explicitly set no algorithm
  default: ReportInvalidAlgorithm(alg);
  }

  if (verboseLevel) {
    G4cout << " >>> G4HadDecayGenerator";
    if (theAlgorithm) G4cout << " using " << theAlgorithm->GetName();
    G4cout << G4endl;
  }
}

G4HadDecayGenerator::G4HadDecayGenerator(G4VHadDecayAlgorithm* alg,
					 G4int verbose)
  : verboseLevel(verbose), theAlgorithm(alg) {
  if (verboseLevel) {
    G4cout << " >>> G4HadDecayGenerator";
    if (theAlgorithm) G4cout << " using " << theAlgorithm->GetName();
    G4cout << G4endl;
  }
}

G4HadDecayGenerator::~G4HadDecayGenerator() {
  delete theAlgorithm;
  theAlgorithm = 0;
}


// Sanity checks -- throws exception if no algorithm chosen

void G4HadDecayGenerator::ReportInvalidAlgorithm(Algorithm alg) const {
  if (verboseLevel) 
    G4cerr << "G4HadDecayGenerator: bad algorithm code " << alg << G4endl;

  throw G4HadronicException(__FILE__, __LINE__, "Invalid algorithm code");
}

void G4HadDecayGenerator::ReportMissingAlgorithm() const {
  if (verboseLevel) 
    G4cerr << "G4HadDecayGenerator: no algorithm specified" << G4endl;

  throw G4HadronicException(__FILE__, __LINE__, "Null algorithm pointer");
}


// Enable (or disable if 0) diagnostic messages
void G4HadDecayGenerator::SetVerboseLevel(G4int verbose) {
  verboseLevel = verbose;
  if (theAlgorithm) theAlgorithm->SetVerboseLevel(verbose);
}

const G4String& G4HadDecayGenerator::GetAlgorithmName() const {
  static const G4String& none = "NONE";
  return (theAlgorithm ? theAlgorithm->GetName() : none);
}


// Initial state (rest mass) and list of final masses

G4bool 
G4HadDecayGenerator::Generate(G4double initialMass,
				   const std::vector<G4double>& masses,
				   std::vector<G4LorentzVector>& finalState) {
  if (verboseLevel) 
    G4cout << " >>> G4HadDecayGenerator::Generate (mass)" << G4endl;

  if (!theAlgorithm) ReportMissingAlgorithm();

  if (masses.size() == 1U)
    return GenerateOneBody(initialMass, masses, finalState);

  theAlgorithm->Generate(initialMass, masses, finalState);
  return !finalState.empty();		// Generator failure returns empty state
}

// Initial state particle and list of final masses

G4bool 
G4HadDecayGenerator::Generate(const G4ParticleDefinition* initialPD,
				   const std::vector<G4double>& masses,
				   std::vector<G4LorentzVector>& finalState) {
  if (verboseLevel) 
    G4cout << " >>> G4HadDecayGenerator::Generate (particle)" << G4endl;

  return (initialPD && Generate(initialPD->GetPDGMass(), masses, finalState));
}

// Final state particles will be boosted to initial-state frame

G4bool 
G4HadDecayGenerator::Generate(const G4LorentzVector& initialState,
				   const std::vector<G4double>& masses,
			      std::vector<G4LorentzVector>& finalState) {
  if (verboseLevel) 
    G4cout << " >>> G4HadDecayGenerator::Generate (frame)" << G4endl;

  G4bool good = Generate(initialState.m(), masses, finalState);
  if (good) {
    G4ThreeVector bv = initialState.boostVector();
    for (size_t i=0; i<finalState.size(); i++) {
      finalState[i].boost(bv);
    }
  }

  return good;
}


// Handle special case of "one body decay" (used for kaon mixing)

G4bool G4HadDecayGenerator::
GenerateOneBody(G4double initialMass,
		const std::vector<G4double>& masses,
		std::vector<G4LorentzVector>& finalState) const {
  if (verboseLevel>1) 
    G4cout << " >>> G4HadDecayGenerator::GenerateOneBody" << G4endl;

  // Initialization and sanity checks
  finalState.clear();

  if (masses.size() != 1U) return false;	// Should not have been called
  if (std::fabs(initialMass-masses[0]) > eV) return false;

  if (verboseLevel>2) G4cout << " finalState mass = " << masses[0] << G4endl;

  finalState.push_back(G4LorentzVector(0.,0.,0.,masses[0]));
  return true;
}
