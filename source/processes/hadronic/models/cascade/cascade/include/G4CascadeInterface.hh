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
// $Id: G4CascadeInterface.hh,v 1.23 2010-09-23 18:13:32 mkelsey Exp $
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

#ifndef G4CASCADEINTERFACE_H
#define G4CASCADEINTERFACE_H 1

#include "G4VIntraNuclearTransportModel.hh"
#include "G4FragmentVector.hh"
#include "G4KineticTrackVector.hh"
#include "G4LorentzRotation.hh"
#include "G4Nucleon.hh"
#include "G4Nucleus.hh"
#include "G4ParticleChange.hh"
#include "G4ReactionProduct.hh"
#include "G4ReactionProductVector.hh"

class G4InuclCollider;
class G4InuclParticle;
class G4CollisionOutput;
class G4CascadeCheckBalance;
class G4V3DNucleus;


class G4CascadeInterface : public G4VIntraNuclearTransportModel {

public:
  G4CascadeInterface(const G4String& name = "Bertini Cascade");

  virtual ~G4CascadeInterface();

  G4ReactionProductVector* Propagate(G4KineticTrackVector* theSecondaries,
				     G4V3DNucleus* theNucleus);
  
  G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
				 G4Nucleus& theNucleus); 

  void setVerboseLevel(G4int verbose) { verboseLevel = verbose; }

  // Select betweeen different post-cascade de-excitation models
  void useCascadeDeexcitation();
  void usePreCompoundDeexcitation();

protected:
  // Convert input projectile and target to Bertini internal types
  void createBullet(const G4HadProjectile& aTrack);
  void createTarget(G4Nucleus& theNucleus);

  // Evaluate whether any outgoing particles penetrated Coulomb barrier
  G4bool coulombBarrierViolation() const;

  // Conditions for rejecting cascade attempt
  G4bool retryInelasticProton() const;
  G4bool retryInelasticNucleus() const;

  // Fill sparse array with minimum momenta for inelastic on hydrogen
  void initializeElasticCuts();

  // Transfer Bertini internal final state to hadronics interface
  void copyOutputToHadronicResult();

  // Terminate job because of energy/momentum/etc. violations
  void throwNonConservationFailure();

private:
  G4int operator==(const G4CascadeInterface& right) const {
    return (this == &right);
  }

  G4int operator!=(const G4CascadeInterface& right) const {
    return (this != &right);
  }

  static const G4int maximumTries;	// Number of iterations for inelastic

  G4double cutElastic[32];		// Bullet momenta for hydrogen target

  G4int verboseLevel;
  G4int numberOfTries;

  G4HadFinalState theResult;
  G4InuclCollider* collider;
  G4CascadeCheckBalance* balance;

  G4InuclParticle* bullet;
  G4InuclParticle* target;
  G4CollisionOutput* output;

  G4LorentzRotation bulletInLabFrame;
};

#endif // G4CASCADEINTERFACE_H
