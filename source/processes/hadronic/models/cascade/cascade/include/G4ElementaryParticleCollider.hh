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
// $Id: G4ElementaryParticleCollider.hh 71940 2013-06-28 19:04:44Z mkelsey $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100315  M. Kelsey -- Remove "using" directive and unnecessary #includes.
// 20100407  M. Kelsey -- Eliminate return-by-value std::vector<> by creating
//		three data buffers for particles, momenta, and particle types.
//		The following functions now return void and are non-const:
//		  ::generateSCMfinalState()
//		  ::generateMomModules() (also remove input vector)
//		  ::generateStrangeChannelPartTypes()
//		  ::generateSCMpionAbsorption()
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide(); merge
//		public vs. private ::collide() functions.
// 20100511  M. Kelsey -- Remove G4PionSampler and G4NucleonSampler.  Expand
//		particle-types selector to all modes, not just strangeness.
// 20100517  M. Kelsey -- Inherit from common base class, make arrays static
// 20100714  M. Kelsey -- Switch to new G4CascadeColliderBase class
// 20100726  M. Kelsey -- Move remaining std::vector<> buffers here
// 20100804  M. Kelsey -- Add printFinalStateTables() function.
// 20110923  M. Kelsey -- Add optional stream& to printFinalStateTables().
// 20130129  M. Kelsey -- Add static arrays and interpolators for two-body
//		angular distributions (addresses MT thread-local issue)
// 20130131  D. Wright -- Use new *AngDst classes for gamma-N two-body
// 20130220  M. Kelsey -- Replace two-body angular code with new *AngDst classes
// 20130221  M. Kelsey -- Move two-body angular dist classes to factory
// 20130306  M. Kelsey -- Move printFinalStateTables() to table factory
// 20130307  M. Kelsey -- Reverse order of dimensions for rmn array
// 20130422  M. Kelsey -- Move kinematics to G4CascadeFinalStateAlgorithm
// 20130508  D. Wright -- Add muon capture, with absorption on quasideuterons
// 20130620  Address Coverity complaint about missing copy actions
// 20130628  Add function to split dibaryon into particle_kinds list
// 20141009  M. Kelsey -- Add pion absorption by single nucleons, with
//		nuclear recoil.  Improves pi- capture performance.

#ifndef G4ELEMENTARY_PARTICLE_COLLIDER_HH
#define G4ELEMENTARY_PARTICLE_COLLIDER_HH

#include "G4CascadeColliderBase.hh"
#include "G4CascadeFinalStateGenerator.hh"
#include "G4CascadeInterpolator.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4LorentzVector.hh"
#include <iosfwd>
#include <vector>

class G4CollisionOutput;


class G4ElementaryParticleCollider : public G4CascadeColliderBase {
public:
  G4ElementaryParticleCollider();
  virtual ~G4ElementaryParticleCollider() {};
  
  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& output);

  void setNucleusState(G4int a, G4int z) {	// Nucleus to use for recoil
    nucleusA = a; nucleusZ = z;
  }

private:
  G4int generateMultiplicity(G4int is, G4double ekin) const;

  void generateOutgoingPartTypes(G4int is, G4int mult, G4double ekin);

  void generateSCMfinalState(G4double ekin, G4double etot_scm,
			     G4InuclElementaryParticle* particle1,
			     G4InuclElementaryParticle* particle2); 

  // Pion (or photon) absorption on a dibaryon
  void generateSCMpionAbsorption(G4double etot_scm,
				 G4InuclElementaryParticle* particle1,
				 G4InuclElementaryParticle* particle2); 

  // Muon absorption on a dibaryon (with outgoing neutrino)
  void generateSCMmuonAbsorption(G4double etot_scm,
				 G4InuclElementaryParticle* particle1,
				 G4InuclElementaryParticle* particle2); 

  // Pion absorption on a single nucleon (charge exchange)
  void generateSCMpionNAbsorption(G4double etot_scm,
				  G4InuclElementaryParticle* particle1,
				  G4InuclElementaryParticle* particle2); 

  G4bool pionNucleonAbsorption(G4double ekin) const;

  G4bool splitQuasiDeuteron(G4int qdtype); 	// Fill kinds with NN components

  void fillOutgoingMasses();		// Fill mass arrays from particle types

  // Utility class to generate final-state kinematics
  G4CascadeFinalStateGenerator fsGenerator;

  // Internal buffers for lists of secondaries
  std::vector<G4InuclElementaryParticle> particles;
  std::vector<G4LorentzVector> scm_momentums;
  std::vector<G4double> modules;
  std::vector<G4double> masses;
  std::vector<G4double> masses2;
  std::vector<G4int> particle_kinds;

  // Nuclear environment (to do pion-nucleon absorption)
  G4int nucleusA, nucleusZ;

private:
  // Copying of modules is forbidden
  G4ElementaryParticleCollider(const G4ElementaryParticleCollider&);
  G4ElementaryParticleCollider& operator=(const G4ElementaryParticleCollider&);
};

#endif	/* G4ELEMENTARY_PARTICLE_COLLIDER_HH */


