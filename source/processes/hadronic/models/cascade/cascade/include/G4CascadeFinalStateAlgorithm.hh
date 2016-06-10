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
// Author:  Michael Kelsey (SLAC)
// Date:    15 April 2013
//
// Description: Subclass of models/util G4VHadDecayAlgorithm which uses
//		old INUCL parametrizations for momentum and angular
//		distributions.
//

#ifndef G4CascadeFinalStateAlgorithm_hh
#define G4CascadeFinalStateAlgorithm_hh 1

#include "globals.hh"
#include "G4VHadDecayAlgorithm.hh"
#include "G4LorentzConvertor.hh"

class G4InuclElementaryParticle;
class G4MultiBodyMomentumDist;
class G4TwoBodyAngularDist;
class G4VMultiBodyMomDst;
class G4VTwoBodyAngDst;


class G4CascadeFinalStateAlgorithm : public G4VHadDecayAlgorithm {
public:
  G4CascadeFinalStateAlgorithm();
  virtual ~G4CascadeFinalStateAlgorithm();

  virtual void SetVerboseLevel(G4int verbose);	// Pass through to factories

  // Select appropriate distributions based on interaction
  void Configure(G4InuclElementaryParticle* bullet,
		 G4InuclElementaryParticle* target,
		 const std::vector<G4int>& particle_kinds);

protected:
  // Two-body generation uses angular-distribution function
  virtual void GenerateTwoBody(G4double initialMass,
			       const std::vector<G4double>& masses,
			       std::vector<G4LorentzVector>& finalState);

  // N-body generation uses momentum-modulus distribution, computed angles
  virtual void GenerateMultiBody(G4double initialMass,
				 const std::vector<G4double>& masses,
				 std::vector<G4LorentzVector>& finalState);

  // Compute kinematic quantities needed for distributions
  void SaveKinematics(G4InuclElementaryParticle* bullet,
		      G4InuclElementaryParticle* target);

  // Select generator based on initial and final state
  void ChooseGenerators(G4int is, G4int fs);

  // Generate momentum magnitudes and validate for use
  void FillMagnitudes(G4double initialMass,
		      const std::vector<G4double>& masses);

  G4bool satisfyTriangle(const std::vector<G4double>& pmod) const;

  // Generate momentum directions into final state
  void FillDirections(G4double initialMass,
		      const std::vector<G4double>& masses,
		      std::vector<G4LorentzVector>& finalState);

  void FillDirThreeBody(G4double initialMass,
			const std::vector<G4double>& masses,
			std::vector<G4LorentzVector>& finalState);

  void FillDirManyBody(G4double initialMass,
		       const std::vector<G4double>& masses,
		       std::vector<G4LorentzVector>& finalState);

  G4double GenerateCosTheta(G4int ptype, G4double pmod) const;

  // SPECIAL:  Generate N-body phase space using Kopylov algorithm
  void FillUsingKopylov(G4double initialMass,
			const std::vector<G4double>& masses,
			std::vector<G4LorentzVector>& finalState);

  G4double BetaKopylov(G4int K) const;	// Copied from G4HadPhaseSpaceKopylov

private:
  const G4VMultiBodyMomDst* momDist;	// Buffers for selected distributions
  const G4VTwoBodyAngDst* angDist;	// Will be NULL for 3+body channels

  std::vector<G4int> kinds;		// Copy of particle_kinds list
  G4int multiplicity;			// Final state size, for convenience
  G4double bullet_ekin;			// Kinematics needed for distributions
  G4LorentzConvertor toSCM;		// Handles complex rotations/transforms

  std::vector<G4double> modules;	// Buffers for generating momenta
  G4ThreeVector mom;

  static const G4double maxCosTheta;	// Cut for valid polar angle generation
  static const G4double oneOverE;	// Numeric value of 1/e for calculations
  static const G4double small;		// Cut for momentum/kinematics
  static const G4int itry_max;		// Maximum number of generation attempts
};

#endif	/* G4CascadeFinalStateAlgorithm_hh */
