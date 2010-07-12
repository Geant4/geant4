#ifndef G4CASCADE_CHECK_BALANCE_HH
#define G4CASCADE_CHECK_BALANCE_HH
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
// $Id: G4CascadeCheckBalance.hh,v 1.4 2010-07-12 05:28:33 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// Verify and report four-momentum conservation for collision output; uses
// same interface as collision generators.
//
// 20100624  M. Kelsey -- Add baryon conservation check and kinetic energy
// 20100628  M. Kelsey -- Add interface to take list of particles directly
// 20100711  M. Kelsey -- Add name of parent Collider for reporting messages,
//		allow changing parent name, add interface for nuclear fragments

#include "G4VCascadeCollider.hh"
#include "globals.hh"
#include "G4LorentzVector.hh"
#include "G4InuclElementaryParticle.hh"
#include <vector>

class G4InuclParticle;
class G4CollisionOutput;

class G4CascadeCheckBalance : public G4VCascadeCollider {
public:
  G4CascadeCheckBalance(G4double relative, G4double absolute,
			const char* owner="G4CascadeCheckBalance");
  virtual ~G4CascadeCheckBalance() {};

  void setOwner(const char* owner) { setName(owner); }

  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& output);

  // This is for use with G4EPCollider internal checks
  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       const std::vector<G4InuclElementaryParticle>& particles);

  // This is for use with G4Fissioner internal checks
  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       const std::vector<G4InuclNuclei>& fragments);

  // Checks on conservation laws (kinematics, baryon number, charge)
  G4bool energyOkay() const;
  G4bool ekinOkay() const;
  G4bool momentumOkay() const;
  G4bool baryonOkay() const;
  G4bool chargeOkay() const;

  // Global check, used by G4CascadeInterface validation loop
  G4bool okay() const { return (energyOkay() && momentumOkay() &&
				baryonOkay() && chargeOkay()); }

  // Calculations of conserved quantities from initial and final state
  G4double deltaE() const { return (final.e() - initial.e()); }
  G4double relativeE() const { return (deltaE()==0.)?0.:deltaE()/initial.e(); }

  G4double deltaKE() const { return (ekin(final) - ekin(initial)); }
  G4double relativeKE() const {
    return (ekin(initial)<=0.) ? 0. : deltaKE()/ekin(initial);
  }

  G4double deltaP() const { return (final.rho() - initial.rho()); }
  G4double relativeP() const { return (deltaP()==0.)?0.:deltaP()/initial.rho(); }

  // Baryon number and charge are discrete; no bounds and no "relative" scale
  G4double deltaB() const { return (finalBaryon - initialBaryon); }
  G4double deltaQ() const { return (finalCharge - initialCharge); }

protected:
  // Utility function for kinetic energy
  G4double ekin(const G4LorentzVector& p) const { return (p.e() - p.m()); }

private:
  G4double relativeLimit;
  G4double absoluteLimit;

  G4LorentzVector initial;	// Four-vectors for computing violations
  G4LorentzVector final;

  G4int initialBaryon;		// Total baryon number
  G4int finalBaryon;

  G4int initialCharge;		// Total charge
  G4int finalCharge;
};

#endif	/* G4CASCADE_CHECK_BALANCE_HH */
