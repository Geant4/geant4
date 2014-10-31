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
#ifndef G4CASCADE_CHECK_BALANCE_HH
#define G4CASCADE_CHECK_BALANCE_HH
// $Id: G4CascadeCheckBalance.hh 71942 2013-06-28 19:08:11Z mkelsey $
//
// Verify and report four-momentum conservation for collision output; uses
// same interface as collision generators.
//
// 20100624  M. Kelsey -- Add baryon conservation check and kinetic energy
// 20100628  M. Kelsey -- Add interface to take list of particles directly
// 20100711  M. Kelsey -- Add name of parent Collider for reporting messages,
//		allow changing parent name, add interface for nuclear fragments
// 20100713  M. Kelsey -- Add (publicly adjustable) tolerance for zeroes
// 20100715  M. Kelsey -- FPE!  Need to check initial values before computing
//		relative error.
// 20100715  M. Kelsey -- Add G4CascadParticle interface for G4NucleiModel;
//		do momentum check on direction, not just magnitude.  Move
//		temporary G4CollisionOutput buffer here, for thread-safety
// 20100909  M. Kelsey -- Add interface to get four-vector difference, and
//		to supply both kinds of particle lists (G4IntraNucleiCascader)
// 20100923  M. Kelsey -- Baryon and charge deltas should have been integer
// 20110328  M. Kelsey -- Add default ctor and explicit limit setting
// 20110722  M. Kelsey -- For IntraNucleiCascader, take G4CollOut as argument
// 20121002  M. Kelsey -- Add strangeness check (useful for Omega- beam)
// 20130620  Address Coverity complaint about missing copy actions
// 20130621  Add interface to take G4Fragment input instead of G4InuclNuclei.
// 20140930  Change name from "const char*" to "const G4String"

#include "G4VCascadeCollider.hh"
#include "globals.hh"
#include "G4CollisionOutput.hh"
#include "G4LorentzVector.hh"
#include <cmath>
#include <vector>

class G4CascadParticle;
class G4InuclElementaryParticle;
class G4InuclNuclei;
class G4InuclParticle;

class G4CascadeCheckBalance : public G4VCascadeCollider {
public:
  static const G4double tolerance;	// Don't do floating zero!

  explicit G4CascadeCheckBalance(const G4String& owner="G4CascadeCheckBalance");

  G4CascadeCheckBalance(G4double relative, G4double absolute,
			const G4String& owner="G4CascadeCheckBalance");
  virtual ~G4CascadeCheckBalance() {};

  void setOwner(const G4String& owner) { setName(owner); }

  void setLimits(G4double relative, G4double absolute) {
    setRelativeLimit(relative);
    setAbsoluteLimit(absolute);
  }

  void setRelativeLimit(G4double limit) { relativeLimit = limit; }
  void setAbsoluteLimit(G4double limit) { absoluteLimit = limit; }

  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& output);

  // This is for use with G4VCascadeDeexcitation modules
  void collide(const G4Fragment& fragment, G4CollisionOutput& output);

  // This is for use with G4EPCollider internal checks
  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       const std::vector<G4InuclElementaryParticle>& particles);

  // This is for use with G4NucleiModel internal checks
  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       const std::vector<G4CascadParticle>& particles);

  // This is for use with G4IntraNucleiCascader
  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& output,
	       const std::vector<G4CascadParticle>& cparticles);

  // This is for use with G4BigBanger internal checks
  void collide(const G4Fragment& target,
	       const std::vector<G4InuclElementaryParticle>& particles);

  // This is for use with G4Fissioner internal checks
  void collide(const G4Fragment& target,
	       const std::vector<G4InuclNuclei>& fragments);

  // Checks on conservation laws (kinematics, baryon number, charge, hyperons)
  G4bool energyOkay() const;
  G4bool ekinOkay() const;
  G4bool momentumOkay() const;
  G4bool baryonOkay() const;
  G4bool chargeOkay() const;
  G4bool strangeOkay() const;

  // Global check, used by G4CascadeInterface validation loop
  // NOTE:  Strangeness is not required to be conserved in final state
  G4bool okay() const { return (energyOkay() && momentumOkay() &&
				baryonOkay() && chargeOkay()); }

  // Calculations of conserved quantities from initial and final state
  // FIXME:  Relative comparisons don't work for zero!
  G4double deltaE() const { return (final.e() - initial.e()); }
  G4double relativeE() const {
    return ( (std::abs(deltaE())<tolerance) ? 0. : 
	     (initial.e()<tolerance) ? 1. : deltaE()/initial.e() );
  }

  G4double deltaKE() const { return (ekin(final) - ekin(initial)); }
  G4double relativeKE() const {
    return ( (std::abs(deltaKE())<tolerance) ? 0. : 
	     (ekin(initial)<tolerance) ? 1. : deltaKE()/ekin(initial) );
  }

  G4double deltaP() const { return deltaLV().rho(); }
  G4double relativeP() const {
    return ( (std::abs(deltaP())<tolerance) ? 0. : 
	     (initial.rho()<tolerance) ? 1. : deltaP()/initial.rho() );
  }

  G4LorentzVector deltaLV() const { return final - initial; }

  // Baryon number, charge, S are discrete; no bounds and no "relative" scale
  G4int deltaB() const { return (finalBaryon - initialBaryon); }
  G4int deltaQ() const { return (finalCharge - initialCharge); }
  G4int deltaS() const { return (finalStrange- initialStrange); }

protected:
  // Utility function for kinetic energy
  G4double ekin(const G4LorentzVector& p) const { return (p.e() - p.m()); }

private:
  G4double relativeLimit;	// Fractional bound on conservation
  G4double absoluteLimit;	// Absolute (GeV) bound on conservation

  G4LorentzVector initial;	// Four-vectors for computing violations
  G4LorentzVector final;

  G4int initialBaryon;		// Total baryon number
  G4int finalBaryon;

  G4int initialCharge;		// Total charge
  G4int finalCharge;

  G4int initialStrange;		// Total strangeness (s-quark content)
  G4int finalStrange;

  G4CollisionOutput tempOutput;		// Buffer for direct-list interfaces

private:
  // Copying of modules is forbidden
  G4CascadeCheckBalance(const G4CascadeCheckBalance&);
  G4CascadeCheckBalance& operator=(const G4CascadeCheckBalance&);
};

#endif	/* G4CASCADE_CHECK_BALANCE_HH */
