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
// $Id: G4IntraNucleiCascader.hh 71719 2013-06-21 00:01:54Z mkelsey $
//
// 20100315  M. Kelsey -- Remove "using" directory and unnecessary #includes.
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100517  M. Kelsey -- Inherit from common base class, make other colliders
//		simple data members
// 20100617  M. Kelsey -- Make G4NucleiModel a data member, instead of
//		creating and deleting on every cycle.
// 20100623  M. Kelsey -- Undo change from 0617.  G4NucleiModel not reusable.
// 20100714  M. Kelsey -- Switch to new G4CascadeColliderBase class
// 20100716  M. Kelsey -- Eliminate inter_case; use base-class functionality,
//		add function to compute recoil nuclear mass on the fly
// 20100720  M. Kelsey -- Make EPCollider pointer member
// 20100722  M. Kelsey -- Move cascade output buffers to .hh file
// 20100728  M. Kelsey -- Move G4NucleiModel here, as pointer member
// 20100907  M. Kelsey -- Add new "makeResidualFragment" to create
//		G4InuclNuclei at current stage of cascade
// 20100909  M. Kelsey -- Drop makeResidualFragment(), getResidualMass() and
//		local G4InuclNuclei object, replace with new RecoilMaker.
//		Move goodCase() to RecoilMaker.
// 20100916  M. Kelsey -- Add functions to handle trapped particles, and to
//		decay hyperons.
// 20110224  M. Kelsey -- Add ::rescatter() function which takes a list of
//		pre-existing secondaries as input.  Split ::collide() into
//		separate utility functions.  Move cascade parameters to static
//		data members.  Add setVerboseLevel().
// 20110303  M. Kelsey -- Add more cascade functions to support rescattering
// 20110304  M. Kelsey -- Modify rescatter to use original Propagate() input
// 20110316  M. Kelsey -- Add function to do G4KineticTrack conversion, decay
//		rescattering resonances in situ.
// 20110324  M. Kelsey -- Add list of nucleon hit locations for rescatter().
// 20110721  M. Kelsey -- Drop decayTrappedParticle(G4KineticTrack*).
// 20110722  M. Kelsey -- Deprecate "output_particles" list in favor of using
//		output directly (will help with pre-cascade issues).
// 20110729  M. Kelsey -- Replace convertKineticToCascade() to reduce churn.
// 20110801  M. Kelsey -- Add local target buffers for rescattering, to avoid
//		memory leak.
// 20110919  M. Kelsey -- Add optional final-state clustering
// 20130304  M. Kelsey -- Add new G4CascadeHistory for cacasde structure reporting
// 20130620  Address Coverity complaint about missing copy actions
// 20141204  M. Kelsey -- Add function to test for non-interacting particles

#ifndef G4INTRA_NUCLEI_CASCADER_HH
#define G4INTRA_NUCLEI_CASCADER_HH

#include "G4CascadeColliderBase.hh"
#include "G4CollisionOutput.hh"
#include "G4ThreeVector.hh"
#include <vector>

class G4CascadParticle;
class G4CascadeCoalescence;
class G4CascadeHistory;
class G4CascadeRecoilMaker;
class G4ElementaryParticleCollider;
class G4InuclElementaryParticle;
class G4InuclParticle;
class G4KineticTrack;
class G4KineticTrackVector;
class G4NucleiModel;
class G4V3DNucleus;


class G4IntraNucleiCascader : public G4CascadeColliderBase {
public:
  G4IntraNucleiCascader();
  virtual ~G4IntraNucleiCascader();

  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& globalOutput);

  // For use with Propagate to preload a set of secondaries
  void rescatter(G4InuclParticle* bullet, G4KineticTrackVector* theSecondaries,
		 G4V3DNucleus* theNucleus, G4CollisionOutput& globalOutput);

  void setVerboseLevel(G4int verbose=0);

private:
  static const G4int itry_max;		// Maximum number of attempts
  static const G4int reflection_cut;	// Maximum number of reflections
  static const G4double small_ekin;	// Tolerance for round-off zero
  static const G4double quasielast_cut;	// To recover elastic scatters

protected:
  G4bool initialize(G4InuclParticle* bullet, G4InuclParticle* target);

  void newCascade(G4int itry);		// Clear buffers for next attempt
  void setupCascade();			// Fill cascade using nuclear model
  void generateCascade();		// Track secondaries through nucleus
  G4bool finishCascade();		// Clean up output, check consistency

  void finalize(G4int itry, 		// Transfer final state for return
		G4InuclParticle* bullet, G4InuclParticle* target,
		G4CollisionOutput& globalOutput);

  G4InuclParticle* createTarget(G4V3DNucleus* theNucleus);

  // Functions to transfer input high-energy cascade for propagation
  void preloadCascade(G4V3DNucleus* theNucleus,
		      G4KineticTrackVector* theSecondaries);
  void copyWoundedNucleus(G4V3DNucleus* theNucleus);
  void copySecondaries(G4KineticTrackVector* theSecondaries);
  void processSecondary(const G4KineticTrack* aSecondary);
  void releaseSecondary(const G4KineticTrack* aSecondary);

  // Functions to handle, e.g., low-energy hyperons stuck inside potential
  void processTrappedParticle(const G4CascadParticle& trapped);
  void decayTrappedParticle(const G4CascadParticle& trapped);

  // Test if particle is able to interact in nucleus
  G4bool particleCanInteract(const G4CascadParticle& cpart) const;

private: 
  G4NucleiModel* model;
  G4ElementaryParticleCollider* theElementaryParticleCollider;
  G4CascadeRecoilMaker* theRecoilMaker;
  G4CascadeCoalescence* theClusterMaker;
  G4CascadeHistory* theCascadeHistory;

  // Buffers and parameters for cascade attempts
  G4InuclNuclei* tnuclei;		// Target nucleus (must be non-zero)
  G4InuclNuclei* bnuclei;		// Non-zero if ion-ion collision
  G4InuclElementaryParticle* bparticle;	// Non-zero if hadron-ion collision

  G4double minimum_recoil_A;		// Require fragment with this mass
  G4double coulombBarrier;

  // Buffers for creation (and reuse) of rescattering targets
  G4InuclNuclei* nucleusTarget;
  G4InuclElementaryParticle* protonTarget;

  // Buffers for collecting result of cascade (reset on each iteration)
  G4CollisionOutput output;
  std::vector<G4CascadParticle> cascad_particles;
  std::vector<G4CascadParticle> new_cascad_particles;
  G4ExitonConfiguration theExitonConfiguration;

  std::vector<G4ThreeVector> hitNucleons;	// Nucleons hit before rescatter

private:
  // Copying of modules is forbidden
  G4IntraNucleiCascader(const G4IntraNucleiCascader&);
  G4IntraNucleiCascader& operator=(const G4IntraNucleiCascader&);
};        

#endif /* G4INTRA_NUCLEI_CASCADER_HH */
