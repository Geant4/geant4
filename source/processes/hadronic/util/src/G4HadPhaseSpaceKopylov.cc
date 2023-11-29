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
// Multibody "phase space" generator using Kopylov's algorithm
//
// Author:	Michael Kelsey (SLAC) <kelsey@slac.stanford.edu>

#include "G4HadPhaseSpaceKopylov.hh"
#include "G4LorentzVector.hh"
#include "G4Pow.hh"
#include "Randomize.hh"
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>


// Generator

void G4HadPhaseSpaceKopylov::
GenerateMultiBody(G4double initialMass,
		  const std::vector<G4double>& masses,
		  std::vector<G4LorentzVector>& finalState) {
  if (GetVerboseLevel()) G4cout << GetName() << "::GenerateMultiBody" << G4endl;

  finalState.clear();

  G4int N = (G4int)masses.size();
  finalState.resize(N);

  G4double mtot = std::accumulate(masses.begin(), masses.end(), 0.0);
  G4double mu = mtot;
  G4double PFragMagCM = 0.0;
  G4double Mass = initialMass;
  G4double T = Mass-mtot;
  G4LorentzVector PFragCM(0.0,0.0,0.0,0.0);
  G4LorentzVector PRestCM(0.0,0.0,0.0,0.0);
  G4LorentzVector PRestLab(0.0,0.0,0.0,Mass);

  for (G4int k=N-1; k>0; --k) {
    mu -= masses[k];
    T *= (k>1) ? BetaKopylov(k) : 0.;
    
    G4double RestMass = mu + T;
    
    PFragMagCM = TwoBodyMomentum(Mass,masses[k],RestMass);
    
    // Create a unit vector with a random direction isotropically distributed
    G4ThreeVector RandVector = UniformVector(PFragMagCM);
    
    PFragCM.setVectM(RandVector,masses[k]);
    PRestCM.setVectM(-RandVector,RestMass);

    G4ThreeVector BoostV = PRestLab.boostVector();
    
    PFragCM.boost(BoostV);
    PRestCM.boost(BoostV);
    PRestLab = PRestCM;
    Mass = RestMass;
    finalState[k] = PFragCM;
  }
  
  finalState[0] = PRestLab;
}


// Generate scale factor for final state particle

G4double G4HadPhaseSpaceKopylov::BetaKopylov(G4int K) const {
  G4Pow* g4pow = G4Pow::GetInstance();

  G4int N = 3*K - 5;
  G4double xN = G4double(N);
  G4double Fmax = std::sqrt(g4pow->powN(xN/(xN+1.),N)/(xN+1.)); 

  G4double F, chi;
  const G4int maxNumberOfLoops = 10000;
  G4int loopCounter = 0;
  do {
    chi = G4UniformRand();
    F = std::sqrt(g4pow->powN(chi,N)*(1.-chi));      
  } while ( ( Fmax*G4UniformRand() > F ) && ++loopCounter < maxNumberOfLoops );  /* Loop checking, 02.11.2015, A.Ribon */ 
  if ( loopCounter >= maxNumberOfLoops ) {
    G4ExceptionDescription ed;
    ed << " Failed sampling after maxNumberOfLoops attempts : forced exit" << G4endl;
    G4Exception( " G4HadPhaseSpaceKopylov::BetaKopylov ", "HAD_KOPYLOV_001", JustWarning, ed );
  }

  return chi;
}
