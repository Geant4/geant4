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

#ifndef G4HadPhaseSpaceGenbod_HH
#define G4HadPhaseSpaceGenbod_HH 1

#include "G4VHadPhaseSpaceAlgorithm.hh"


class G4HadPhaseSpaceGenbod : public G4VHadPhaseSpaceAlgorithm {
public:
  G4HadPhaseSpaceGenbod(G4int verbose=0);
  virtual ~G4HadPhaseSpaceGenbod() {;}

protected:
  virtual void GenerateMultiBody(G4double initialMass,
				 const std::vector<G4double>& masses,
				 std::vector<G4LorentzVector>& finalState);

protected:
  void Initialize(G4double initialMass,
		  const std::vector<G4double>& masses);
			  
  void FillRandomBuffer();

  void ComputeWeightScale(const std::vector<G4double>& masses);

  void FillEnergySteps(G4double initialMass,
		       const std::vector<G4double>& masses);

  void GenerateMomenta(const std::vector<G4double>& masses,
		       std::vector<G4LorentzVector>& finalState);

  void AccumulateFinalState(size_t i,
			    const std::vector<G4double>& masses,
			    std::vector<G4LorentzVector>& finalState);

  G4bool AcceptEvent() const;	// Use accept-reject to generate distribution
  G4double ComputeWeight() const;

private:
  size_t nFinal;		// Number of final state particles
  G4double totalMass;		// Sum of final state masses
  G4double massExcess;		// Available kinetic energy
  G4double weightMax;		// Maximum possible weight
  G4int nTrials;		// Accept/reject cycles taken

  std::vector<G4double> msum;	// Cumulative sum of masses
  std::vector<G4double> msq;	// Final state squared masses
  std::vector<G4double> rndm;	// Random sequence for effective masses
  std::vector<G4double> meff;	// Random final-state effective masses
  std::vector<G4double> pd;	// Random momentum magnitudes
};

#endif	/* G4HadPhaseSpaceGenbod_HH */
