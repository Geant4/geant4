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
// Abstract base class for multibody "phase space" generators.  Subclasses
// implement a specific algorithm, such as Kopylov, GENBOD, or Makoto's
// NBody.  Subclasses are used by G4HadPhaseSpaceGenerator.
//
// Author:	Michael Kelsey (SLAC) <kelsey@slac.stanford.edu>

#include "G4VHadDecayAlgorithm.hh"
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


// Initial state (rest mass) and list of final masses

void G4VHadDecayAlgorithm::Generate(G4double initialMass,
				    const std::vector<G4double>& masses,
				    std::vector<G4LorentzVector>& finalState) {
  if (verboseLevel) G4cout << GetName() << "::Generate" << G4endl;

  // Initialization and sanity check
  finalState.clear();
  if (!IsDecayAllowed(initialMass, masses)) return;

  // Allow different procedures for two-body or N-body distributions
  if (masses.size() == 2U) 
    GenerateTwoBody(initialMass, masses, finalState);
  else 
    GenerateMultiBody(initialMass, masses, finalState);
}


// Base class does very simple validation of configuration

G4bool G4VHadDecayAlgorithm::
IsDecayAllowed(G4double initialMass,
	       const std::vector<G4double>& masses) const {
  G4bool okay =
    (initialMass > 0. && masses.size() >= 2 &&
     initialMass >= std::accumulate(masses.begin(),masses.end(),0.));

  if (verboseLevel) {
    G4cout << GetName() << "::IsDecayAllowed? initialMass " << initialMass
	   << " " << masses.size() << " masses sum "
	   << std::accumulate(masses.begin(),masses.end(),0.) << G4endl;

    if (verboseLevel>1) PrintVector(masses," ",G4cout);

    G4cout << " Returning " << okay << G4endl;
  }

  return okay;
}


// Momentum function (c.f. PDK() function from CERNLIB W515)

G4double G4VHadDecayAlgorithm::TwoBodyMomentum(G4double M0, G4double M1,
					       G4double M2) const {
  G4double PSQ = (M0+M1+M2)*(M0+M1-M2)*(M0-M1+M2)*(M0-M1-M2);
  if (PSQ < 0.) {
    G4cout << GetName() << ":  problem of decay of M(GeV) " << M0/GeV 
	   << " to M1(GeV) " << M1/GeV << " and M2(GeV) " << M2/GeV
	   << " PSQ(MeV) " << PSQ/MeV << " < 0" << G4endl;
    // exception only if the problem is numerically significant
    if (PSQ < -CLHEP::eV) {
      throw G4HadronicException(__FILE__, __LINE__,"Error in decay kinematics");
    }

    PSQ = 0.;
  }

  return std::sqrt(PSQ)/(2.*M0);
}

// Convenience functions for uniform angular distributions

G4double G4VHadDecayAlgorithm::UniformTheta() const {
  return std::acos(2.0*G4UniformRand() - 1.0);
}

G4double G4VHadDecayAlgorithm::UniformPhi() const {
  return twopi*G4UniformRand();
}


// Dump contents of vector to output

void G4VHadDecayAlgorithm::
PrintVector(const std::vector<G4double>& v,
	    const G4String& vname, std::ostream& os) const {
  os << " " << vname << "(" << v.size() << ") ";
  std::copy(v.begin(), v.end(), std::ostream_iterator<G4double>(os, " "));
  os << std::endl;
}
