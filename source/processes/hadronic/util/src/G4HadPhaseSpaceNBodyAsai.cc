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
// Multibody "phase space" generator using Makoto Asai's NBody method.
//
// Author:	Michael Kelsey (SLAC) <kelsey@slac.stanford.edu>

#include "G4HadPhaseSpaceNBodyAsai.hh"
#include "G4LorentzVector.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "Randomize.hh"
#include <algorithm>
#include <functional>
#include <iterator>
#include <numeric>
#include <vector>


namespace {
  // This wraps the existing #define in a true function
  G4double uniformRand() { return G4UniformRand(); }
}


void G4HadPhaseSpaceNBodyAsai::
GenerateMultiBody(G4double initialMass,
		  const std::vector<G4double>& masses,
		  std::vector<G4LorentzVector>& finalState) {
  if (GetVerboseLevel()) G4cout << GetName() << "::GenerateMultiBody" << G4endl;

  finalState.clear();

  //daughters' mass
  G4int numberOfDaughters = (G4int)masses.size();
  G4double sumofmasses =
    std::accumulate(masses.begin(), masses.end(), 0.);
  
  //Calculate daughter momentum
  std::vector<G4double> daughtermomentum(numberOfDaughters);
  std::vector<G4double> sm(numberOfDaughters);
  G4double tmas;
  G4double weight = 1.0;
  G4int numberOfTry = 0; 
  G4int i;

  std::vector<G4double> rd(numberOfDaughters);
  do {
    //Generate random number in descending order
    rd[0] = 1.0;
    std::generate(rd.begin()+1, rd.end(), uniformRand);
    std::sort(rd.begin(), rd.end(), std::greater<G4double>());

    if (GetVerboseLevel()>1) PrintVector(rd,"rd",G4cout);

    //calcurate virtual mass 
    tmas = initialMass -  sumofmasses;
    G4double temp = sumofmasses; 
    for(i =0; i < numberOfDaughters; i++) {
      sm[i] = rd[i]*tmas + temp;
      temp -= masses[i];
      if (GetVerboseLevel()>1) {
        G4cout << i << " random number:" << rd[i]
	       << " virtual mass:" << sm[i]/GeV << " GeV/c2" <<G4endl; 
      }
    }

    //Calculate daughter momentum
    weight = 1.0;
    i = numberOfDaughters-1;
    daughtermomentum[i] = TwoBodyMomentum(sm[i-1],masses[i-1],sm[i]);
    if (GetVerboseLevel()>1) {
      G4cout << " daughter " << i << ": momentum "
	     << daughtermomentum[i]/GeV << " GeV/c" <<G4endl;
    }
    for(i =numberOfDaughters-2; i>=0; i--) {
      // calculate 
      daughtermomentum[i] = TwoBodyMomentum(sm[i],masses[i],sm[i+1]);
      if(daughtermomentum[i] < 0.0) {
        // !!! illegal momentum !!!
        if (GetVerboseLevel()>0) {
          G4cout << "G4HadPhaseSpaceNBodyAsai::Generate "
		 << " can not calculate daughter momentum "
		 << "\n initialMass " << initialMass/GeV << " GeV/c2"
		 << "\n daughter " << i << ": mass "
		 << masses[i]/GeV << " GeV/c2; momentum "
		 << daughtermomentum[i]/GeV << " GeV/c" << G4endl;
        }
	return;   		// Error detection
      }

      // calculate weight of this events
      weight *= daughtermomentum[i]/sm[i];
      if (GetVerboseLevel()>1) {
	G4cout << " daughter " << i << ": momentum "
	       << daughtermomentum[i]/GeV << " GeV/c" <<G4endl;
      }
    }
    if (GetVerboseLevel()>1) {
      G4cout << " weight: " << weight <<G4endl;
    }
    
    // exit if number of Try exceeds 100
    if (numberOfTry++ > 100) {
      if (GetVerboseLevel()>0) {
        G4cout << "G4HadPhaseSpaceNBodyAsai::Generate "
	       << " can not determine Decay Kinematics " << G4endl;
      }
      return;			// Error detection
    }
  } while (weight > G4UniformRand());  /* Loop checking, 02.11.2015, A.Ribon */

  if (GetVerboseLevel()>1) {
      G4cout << "Start calculation of daughters momentum vector "<<G4endl;
  }
  
  G4double beta;

  finalState.resize(numberOfDaughters);

  i = numberOfDaughters-2;

  G4ThreeVector direction = UniformVector(daughtermomentum[i]);

  finalState[i].setVectM(direction, masses[i]);
  finalState[i+1].setVectM(-direction, masses[i+1]);

  for (i = numberOfDaughters-3;  i >= 0; i--) {
    direction = UniformVector();

    //create daughter particle
    finalState[i].setVectM(-daughtermomentum[i]*direction, masses[i]);

    // boost already created particles 
    beta = daughtermomentum[i];
    beta /= std::sqrt(beta*beta + sm[i+1]*sm[i+1]);
    for (G4int j = i+1; j<numberOfDaughters; j++) {
      finalState[j].boost(beta*direction);
    }
  }
}
