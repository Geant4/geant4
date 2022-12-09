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
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100407  M. Kelsey -- Replace ::resize(0) with ::clear()
// 20100409  M. Kelsey -- Move function code to .cc files, not inlinable
// 20100418  M. Kelsey -- Add function to boost output lists to lab frame
// 20100520  M. Kelsey -- Add function to rotate Z axis, from G4Casc.Interface
// 20100620  M. Kelsey -- Add setVerboseLevel() function
// 20100715  M. Kelsey -- Add total charge and baryon number functions, and a
//		combined "add()" function to put two of these together.
// 20100716  M. Kelsey -- Add interface to handle G4CascadParticles
// 20100924  M. Kelsey -- Use "OutgoingNuclei" name consistently, replacing
//		old "TargetFragment".  Add new (reusable) G4Fragment buffer 
//		and access functions for initial post-cascade processing.
//		Move implementation of add() to .cc file.
// 20100925  M. Kelsey -- Add function to process G4ReactionProduct list
// 20110225  M. Kelsey -- Add interface to remove entries from lists
// 20110311  M. Kelsey -- Add function to boost individual four-vector
// 20110323  M. Kelsey -- Add non-const access to lists (for G4NucleiModel)
// 20110922  M. Kelsey -- Add optional stream argument to printCollisionOutput
// 20121002  M. Kelsey -- Add strangeness calculation
// 20130628  M. Kelsey -- Support multiple recoil fragments (for G4Fissioner)
// 20141208  M. Kelsey -- Split function to do pair-wise "hard" tuning

#ifndef G4COLLISION_OUTPUT_HH
#define G4COLLISION_OUTPUT_HH

#include "G4Fragment.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4LorentzRotation.hh"
#include "G4ReactionProductVector.hh"
#include "G4ios.hh"
#include <iosfwd>
#include <algorithm>
#include <vector>

class G4CascadParticle;
class G4LorentzConvertor;

class G4CollisionOutput
{
 public:

  G4CollisionOutput();
  G4CollisionOutput& operator=(const G4CollisionOutput& right);

  void setVerboseLevel(G4int verbose) { verboseLevel = verbose; };

  // ===== Accumulate contents of lists =====

  void reset();		// Empties lists for new event

  void add(const G4CollisionOutput& right);	// Merge complete objects

  void addOutgoingParticle(const G4InuclElementaryParticle& particle) {
    outgoingParticles.push_back(particle);
  }

  void addOutgoingParticles(const std::vector<G4InuclElementaryParticle>& particles);

  void addOutgoingNucleus(const G4InuclNuclei& nuclei) {
    outgoingNuclei.push_back(nuclei);
  };

  void addOutgoingNuclei(const std::vector<G4InuclNuclei>& nuclea);

  // These are primarily for G4IntraNucleiCascader internal checks
  void addOutgoingParticle(const G4CascadParticle& cparticle);
  void addOutgoingParticles(const std::vector<G4CascadParticle>& cparticles);

  void addOutgoingParticles(const G4ReactionProductVector* rproducts);

  // Special buffer for initial, possible unstable fragments from cascade
  void addRecoilFragment(const G4Fragment* aFragment) {
    if (aFragment) addRecoilFragment(*aFragment);
  }

  void addRecoilFragment(const G4Fragment& aFragment) {
    recoilFragments.push_back(aFragment);
  }
  
  // ===== Remove contents of lists, by index, reference or value  =====

  void removeOutgoingParticle(G4int index);
  void removeOutgoingParticle(const G4InuclElementaryParticle& particle);
  void removeOutgoingParticle(const G4InuclElementaryParticle* particle) {
    if (particle) removeOutgoingParticle(*particle);
  }

  void removeOutgoingNucleus(G4int index);
  void removeOutgoingNucleus(const G4InuclNuclei& nuclei);
  void removeOutgoingNucleus(const G4InuclNuclei* nuclei) {
    if (nuclei) removeOutgoingNucleus(*nuclei);
  }

  void removeRecoilFragment(G4int index=-1);	// No argument removes all

  // ===== Access contents of lists =====

  G4int numberOfOutgoingParticles() const { return (G4int)outgoingParticles.size(); }
    
  const std::vector<G4InuclElementaryParticle>& getOutgoingParticles() const {
    return outgoingParticles;
  };

  std::vector<G4InuclElementaryParticle>& getOutgoingParticles() {
    return outgoingParticles;
  };

  G4int numberOfOutgoingNuclei() const { return (G4int)outgoingNuclei.size(); };
 
  const std::vector<G4InuclNuclei>& getOutgoingNuclei() const {
    return outgoingNuclei;
  };

  std::vector<G4InuclNuclei>& getOutgoingNuclei() { return outgoingNuclei; };

  G4int numberOfFragments() const { return (G4int)recoilFragments.size(); }

  const G4Fragment& getRecoilFragment(G4int index=0) const;

  const std::vector<G4Fragment>& getRecoilFragments() const {
    return recoilFragments;
  };

  std::vector<G4Fragment>& getRecoilFragments() { return recoilFragments; };

  // ===== Get event totals for conservation checking, recoil, etc. ======

  G4LorentzVector getTotalOutputMomentum() const;
  G4int getTotalCharge() const;			// NOTE:  No fractional charges!
  G4int getTotalBaryonNumber() const;
  G4int getTotalStrangeness() const;

  void printCollisionOutput(std::ostream& os=G4cout) const;

  // ===== Manipulate final-state particles for kinematics =====

  void boostToLabFrame(const G4LorentzConvertor& convertor);
  G4LorentzVector boostToLabFrame(G4LorentzVector mom,	// Note pass by value!
				  const G4LorentzConvertor& convertor) const;

  void rotateEvent(const G4LorentzRotation& rotate);
  void trivialise(G4InuclParticle* bullet, G4InuclParticle* target);
  void setOnShell(G4InuclParticle* bullet, G4InuclParticle* target);
  void setRemainingExitationEnergy();

  G4double getRemainingExitationEnergy() const { return eex_rest; };
  G4bool acceptable() const { return on_shell; };

 private: 

  G4int verboseLevel;

  std::vector<G4InuclElementaryParticle> outgoingParticles;
  std::vector<G4InuclNuclei> outgoingNuclei;
  std::vector<G4Fragment> recoilFragments;
  static const G4Fragment emptyFragment;	// To return if list empty

  std::pair<std::pair<G4int,G4int>, G4int> selectPairToTune(G4double de) const; 
  G4bool tuneSelectedPair(G4LorentzVector& mom1, G4LorentzVector& mom2,
			  G4int mom_index) const;

  G4double eex_rest;		// Used by setOnShell() for kinematics
  G4LorentzVector mom_non_cons;
  G4bool on_shell;
};        

#endif // G4COLLISION_OUTPUT_HH 

