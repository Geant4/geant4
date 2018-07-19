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
// $Id: G4CascadeInterface.hh 71938 2013-06-28 19:01:00Z mkelsey $
// Defines an interface to Bertini (BERT) cascade
// based on INUCL  intra-nuclear transport.models 
// with bullet hadron energy ~< 10 GeV
//
// 20100405  M. Kelsey -- Fix constness of op== and op!=
// 20100519  M. Kelsey -- Remove Collider data members
// 20100617  M. Kelsey -- Make G4InuclCollider a local data member
// 20100723  M. Kelsey -- Move G4CollisionOutput here for reuse
// 20100916  M. Kelsey -- Add functions to encapsulate ApplyYourself() actions,
//		make colliders pointers (don't expose dependencies)
// 20100922  M. Kelsey -- Add functions to select de-excitation method
// 20110224  M. Kelsey -- Add createTarget() for use with Propagate(); split
//		conservation law messages to separate function.  Move verbose
//		setting to .cc file, and apply to all member objects.
// 20110301  M. Kelsey -- Add copyPreviousCascade() for use with Propagate()  
//		along with new buffers and related particle-conversion  
//		functions.  Encapulate buffer deletion in clear()
// 20110303  M. Kelsey -- Change "bulletList" name to "inputFragments"
// 20110304  M. Kelsey -- Drop conversion of Propagate() arguments; pass
//		directly to collider for processing.  Rename makeReactionProduct
//		to makeDynamicParticle.
// 20110502  M. Kelsey -- Add filename string to capture random seeds.
// 20110720  M. Kelsey -- Discard elastic-cut array (no longer needed),
//		discard local "theFinalState" (avail in base class).
// 20110801  M. Kelsey -- Make bullet and target buffers local objects (with
//		hadron and nucleus versions) to reduce memory churn
// 20120522  M. Kelsey -- Implement base class IsApplicable, and add overloaded
//		version which takes G4ParticleDefintion, a la G4VProcess.
// 20120822  M. Kelsey -- Add function to dump user configuration settings.
//		Remove local verboseLevel; shadows base class data member.
// 20130501  M. Kelsey -- Add static initializer to created shared objects.
// 20130628  M. Kelsey -- Address Coverity warnings about copy operations.
// 20140116  M. Kelsey -- Move statics to const data members to avoid weird
//		interactions with MT.

#ifndef G4CASCADEINTERFACE_H
#define G4CASCADEINTERFACE_H 1

#include "G4VIntraNuclearTransportModel.hh"
#include "G4FragmentVector.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4LorentzRotation.hh"
#include "G4Nucleon.hh"
#include "G4Nucleus.hh"
#include "G4ParticleChange.hh"
#include "G4ReactionProduct.hh"
#include "G4ReactionProductVector.hh"
#include <vector>

class G4CascadParticle;
class G4CascadeCheckBalance;
class G4CollisionOutput;
class G4DynamicParticle;
class G4HadFinalState;
class G4InuclCollider;
class G4InuclParticle;
class G4ParticleDefinition;
class G4V3DNucleus;


class G4CascadeInterface : public G4VIntraNuclearTransportModel {

public:
  G4CascadeInterface(const G4String& name = "BertiniCascade");

  virtual ~G4CascadeInterface();

  G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries,
				     G4V3DNucleus* theNucleus);
  
  G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
				 G4Nucleus& theNucleus);

  void SetVerboseLevel(G4int verbose);		// Overrides base class

  G4bool IsApplicable(const G4HadProjectile& aTrack,
		      G4Nucleus& theNucleus);

  G4bool IsApplicable(const G4ParticleDefinition* aPD) const;

  // Used with multithreaded applications to preload shared objects
  static void Initialize();

  // Select betweeen different post-cascade de-excitation models
  void useCascadeDeexcitation();
  void usePreCompoundDeexcitation();

  virtual void ModelDescription(std::ostream& outFile) const;
  virtual void DumpConfiguration(std::ostream& outFile) const;

protected:
  void clear();			// Delete previously created particles

  // Convert input projectile and target to Bertini internal types
  G4bool createBullet(const G4HadProjectile& aTrack);
  G4bool createTarget(G4Nucleus& theNucleus);
  G4bool createTarget(G4V3DNucleus* theNucleus);
  G4bool createTarget(G4int A, G4int Z);

  // Evaluate whether any outgoing particles penetrated Coulomb barrier
  G4bool coulombBarrierViolation() const;

  // Conditions for rejecting cascade attempt
  G4bool retryInelasticProton() const;
  G4bool retryInelasticNucleus() const;

  // Transfer Bertini internal final state to hadronics interface
  void copyOutputToHadronicResult();
  G4ReactionProductVector* copyOutputToReactionProducts();

  // Replicate input particles onto output
  G4HadFinalState* NoInteraction(const G4HadProjectile& aTrack,
				 G4Nucleus& theNucleus);

  // Report violations of conservation laws in original frame
  void checkFinalResult();

  // Terminate job because of energy/momentum/etc. violations
  void throwNonConservationFailure();

  // Convert between Bertini and external particle types
  G4DynamicParticle* makeDynamicParticle(const G4InuclElementaryParticle& iep) const;
  G4DynamicParticle* makeDynamicParticle(const G4InuclNuclei& inuc) const;


private:
  G4int operator==(const G4CascadeInterface& right) const {
    return (this == &right);
  }

  G4int operator!=(const G4CascadeInterface& right) const {
    return (this != &right);
  }

  const G4String randomFile;		// Filename to capture random seeds
  const G4int maximumTries;		// Number of iterations for inelastic

  G4int numberOfTries;

  G4InuclCollider* collider;
  G4CascadeCheckBalance* balance;

  G4InuclParticle* bullet;		// Pointers to last filled versions
  G4InuclParticle* target;

  G4CollisionOutput* output;

  G4InuclElementaryParticle hadronBullet;	// Buffers for bullet, target
  G4InuclNuclei             nucleusBullet;

  G4InuclElementaryParticle hadronTarget;
  G4InuclNuclei             nucleusTarget;


private:
  // Copying of modules is forbidden
  G4CascadeInterface(const G4CascadeInterface&);
  G4CascadeInterface& operator=(const G4CascadeInterface&);
};

#endif // G4CASCADEINTERFACE_H
