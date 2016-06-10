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
// Multibody "phase space" generator using GENBOD (CERNLIB W515) method.
//
// Author:	Michael Kelsey (SLAC) <kelsey@slac.stanford.edu>

#include "G4HadPhaseSpaceGenbod.hh"
#include "G4LorentzVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <algorithm>
#include <functional>
#include <iterator>
#include <numeric>
#include <vector>


namespace {
  // Wrap #define in a true function, for passing to std::fill
  G4double uniformRand() { return G4UniformRand(); }
}


// Constructor initializes everything to zero

G4HadPhaseSpaceGenbod::G4HadPhaseSpaceGenbod(G4int verbose)
  : G4VHadPhaseSpaceAlgorithm("G4HadPhaseSpaceGenbod",verbose),
    nFinal(0), totalMass(0.), massExcess(0.), weightMax(0.), nTrials(0) {;}


// C++ re-implementation of GENBOD.F (Raubold-Lynch method)

void G4HadPhaseSpaceGenbod::
GenerateMultiBody(G4double initialMass,
		  const std::vector<G4double>& masses,
		  std::vector<G4LorentzVector>& finalState) {
  if (GetVerboseLevel()) G4cout << GetName() << "::GenerateMultiBody" << G4endl;

  finalState.clear();

  Initialize(initialMass, masses);

  const G4int maxNumberOfLoops = 10000;
  nTrials = 0;
  do {				// Apply accept/reject to get distribution
    ++nTrials;
    FillRandomBuffer();
    FillEnergySteps(initialMass, masses);
  } while ( (!AcceptEvent()) && nTrials < maxNumberOfLoops );  /* Loop checking, 02.11.2015, A.Ribon */
  if ( nTrials >= maxNumberOfLoops ) {
    G4ExceptionDescription ed;
    ed << " Failed sampling after maxNumberOfLoops attempts : forced exit" << G4endl;
    G4Exception( " G4HadPhaseSpaceGenbod::GenerateMultiBody ", "HAD_GENBOD_001", FatalException, ed );
  }
  GenerateMomenta(masses, finalState);
}

void G4HadPhaseSpaceGenbod::
Initialize(G4double initialMass, const std::vector<G4double>& masses) {
  if (GetVerboseLevel()>1) G4cout << GetName() << "::Initialize" << G4endl;

  nFinal = masses.size();
  msum.resize(nFinal, 0.);		// Initialize buffers for filling
  msq.resize(nFinal, 0.);

  std::partial_sum(masses.begin(), masses.end(), msum.begin());
  std::transform(masses.begin(), masses.end(), masses.begin(), msq.begin(),
		 std::multiplies<G4double>());
  totalMass  = msum.back();
  massExcess = initialMass - totalMass;

  if (GetVerboseLevel()>2) {
    PrintVector(msum, "msum", G4cout);
    PrintVector(msq, "msq", G4cout);
    G4cout << " totalMass " << totalMass << " massExcess " << massExcess
	   << G4endl;
  }

  ComputeWeightScale(masses);
}


// Generate ordered list of random numbers

void G4HadPhaseSpaceGenbod::FillRandomBuffer() {
  if (GetVerboseLevel()>1) G4cout << GetName() << "::FillRandomBuffer" << G4endl;

  rndm.resize(nFinal-2,0.);	// Final states generated in sorted order
  std::generate(rndm.begin(), rndm.end(), uniformRand);
  std::sort(rndm.begin(), rndm.end());
  if (GetVerboseLevel()>2) PrintVector(rndm, "rndm", G4cout);
}


// Final state effective masses, min to max

void 
G4HadPhaseSpaceGenbod::FillEnergySteps(G4double initialMass,
				       const std::vector<G4double>& masses) {
  if (GetVerboseLevel()>1) G4cout << GetName() << "::FillEnergySteps" << G4endl;

  meff.clear();
  pd.clear();

  meff.push_back(masses[0]);
  for (size_t i=1; i<nFinal-1; i++) {
    meff.push_back(rndm[i-1]*massExcess + msum[i]);
    pd.push_back(TwoBodyMomentum(meff[i], meff[i-1], masses[i]));
  }
  meff.push_back(initialMass);
  pd.push_back(TwoBodyMomentum(meff[nFinal-1], meff[nFinal-2], masses[nFinal-1]));

  if (GetVerboseLevel()>2) {
    PrintVector(meff,"meff",G4cout);
    PrintVector(pd,"pd",G4cout);
  }
}


// Maximum possible weight for final state (used with accept/reject)

void 
G4HadPhaseSpaceGenbod::ComputeWeightScale(const std::vector<G4double>& masses) {
  if (GetVerboseLevel()>1) 
    G4cout << GetName() << "::ComputeWeightScale" << G4endl;

  weightMax = 1.;
  for (size_t i=1; i<nFinal; i++) {
    weightMax *= TwoBodyMomentum(massExcess+msum[i], msum[i-1], masses[i]);
  }

  if (GetVerboseLevel()>2) G4cout << " weightMax = " << weightMax << G4endl;
}


// Event weight computed as either constant or Fermi-dependent cross-section

G4double G4HadPhaseSpaceGenbod::ComputeWeight() const {
  if (GetVerboseLevel()>1) G4cout << GetName() << "::ComputeWeight" << G4endl;

  return (std::accumulate(pd.begin(), pd.end(), 1./weightMax,
			  std::multiplies<G4double>()));
}

G4bool G4HadPhaseSpaceGenbod::AcceptEvent() const {
  if (GetVerboseLevel()>1) 
    G4cout << GetName() << "::AcceptEvent? " << nTrials << G4endl;

  return (G4UniformRand() <= ComputeWeight());
}


// Final state momentum vectors in CMS system, using Raubold-Lynch method

void G4HadPhaseSpaceGenbod::
GenerateMomenta(const std::vector<G4double>& masses,
		std::vector<G4LorentzVector>& finalState) {
  if (GetVerboseLevel()>1) G4cout << GetName() << "::GenerateMomenta" << G4endl;

  finalState.resize(nFinal);	// Preallocate vectors for convenience below

  for (size_t i=0; i<nFinal; i++) {
    AccumulateFinalState(i, masses, finalState);
    if (GetVerboseLevel()>2)
      G4cout << " finalState[" << i << "] " << finalState[i] << G4endl;
  }
}

// Process final state daughters up to current index

void G4HadPhaseSpaceGenbod::
AccumulateFinalState(size_t i, 
		     const std::vector<G4double>& masses,
		     std::vector<G4LorentzVector>& finalState) {
  if (GetVerboseLevel()>2)
    G4cout << GetName() << "::AccumulateFinalState " << i << G4endl;

  if (i==0) {			// First final state particle left alone
    finalState[i].setVectM(G4ThreeVector(0.,pd[i],0.),masses[i]);
    return;
  }
  
  finalState[i].setVectM(G4ThreeVector(0.,-pd[i-1],0.),masses[i]);
  G4double phi = G4UniformRand() * twopi;
  G4double theta = std::acos(2.*G4UniformRand() - 1.);

  if (GetVerboseLevel() > 2) {
    G4cout << " initialized Py " << -pd[i-1] << " phi " << phi
	   << " theta " << theta << G4endl;
  }

  G4double esys=0.,beta=0.,gamma=1.;
  if (i < nFinal-1) {			// Do not boost final particle
    esys = std::sqrt(pd[i]*pd[i]+meff[i]*meff[i]);
    beta = pd[i] / esys;
    gamma = esys / meff[i];

    if (GetVerboseLevel()>2)
      G4cout << " esys " << esys << " beta " << beta << " gamma " << gamma
	     << G4endl;
  }

  for (size_t j=0; j<=i; j++) {		// Accumulate rotations
    finalState[j].rotateZ(theta).rotateY(phi);
    finalState[j].setY(gamma*(finalState[j].y() + beta*finalState[j].e()));
    if (GetVerboseLevel()>2) G4cout << " j " << j << " " << finalState[j] << G4endl;
  }
}
