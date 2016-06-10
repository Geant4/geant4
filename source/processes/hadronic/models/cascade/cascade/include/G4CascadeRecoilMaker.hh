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
#ifndef G4CASCADE_RECOIL_MAKER_HH
#define G4CASCADE_RECOIL_MAKER_HH
// $Id: G4CascadeRecoilMaker.hh 71719 2013-06-21 00:01:54Z mkelsey $
//
// Collects generated cascade data (using Collider::collide() interface)
// and computes the nuclear recoil kinematics needed to balance the event.
//
// 20100909  M. Kelsey -- Inspired by G4CascadeCheckBalance
// 20100909  M. Kelsey -- Move G4IntraNucleiCascader::goodCase() here, add
//		tolerance for "almost zero" excitation energy, add function
//		to manually override excitation energy
// 20100910  M. Kelsey -- Drop getRecoilFragment() in favor of user calling
//		makeRecoilFragment() with returned non-const pointer.  Drop
//		handling of excitons.
// 20100914  M. Kelsey -- Migrate to integer A and Z
// 20100921  M. Kelsey -- Return G4InuclNuclei using "makeRecoilNuclei()".
//		Repurpose "makeRecoilFragment()" to return G4Fragment.
// 20100924  M. Kelsey -- Add raw excitation energy (mass difference) function
// 20110214  M. Kelsey -- Replace "model" with G4InuclParticle::Model enum
// 20110722  M. Kelsey -- For IntraNucleiCascader, take G4CollOut as argument
// 20130620  Address Coverity complaint about missing copy actions

#include <cmath>
#include <vector>
#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VCascadeCollider.hh"
#include "globals.hh"
#include "G4CollisionOutput.hh"
#include "G4InuclNuclei.hh"
#include "G4Fragment.hh"
#include "G4LorentzVector.hh"

class G4CascadParticle;
class G4CascadeCheckBalance;
class G4InuclElementaryParticle;
class G4InuclParticle;


class G4CascadeRecoilMaker : public G4VCascadeCollider {
public:
  explicit G4CascadeRecoilMaker(G4double tolerance=0.001*CLHEP::MeV);
  virtual ~G4CascadeRecoilMaker();

  // Standard Collider interface (non-const output "buffer")
  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& output);

  // This is for use with G4IntraNucleiCascader
  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& output,
	       const std::vector<G4CascadParticle>& cparticles);

  // Modifiable parameters
  void setTolerance(G4double tolerance) { excTolerance = tolerance; }

  void setRecoilExcitation(G4double Eexc) { excitationEnergy = Eexc; }

  // Build nucleus from current parameters, if physically reasonable
  G4InuclNuclei* makeRecoilNuclei(G4InuclParticle::Model model=G4InuclParticle::DefaultModel);
  G4Fragment* makeRecoilFragment();	// For use with PreCompound

  // Attach exciton configuration for use by "nucleus makers"
  void addExcitonConfiguration(const G4ExitonConfiguration exciton) {
    theExcitons = exciton;
  }

  // Access nuclear configuration parameters
  G4int getRecoilA() const { return recoilA; }
  G4int getRecoilZ() const { return recoilZ; }
  G4double getRecoilExcitation() const { return excitationEnergy; }
  const G4LorentzVector& getRecoilMomentum() const { return recoilMomentum; }

  // Data quality checks
  G4bool goodFragment() const;		// Verify A, Z both meaningful
  G4bool goodRecoil() const;		// And sensible four-vector
  G4bool wholeEvent() const;		// Zero recoil
  G4bool unphysicalRecoil() const { return !wholeEvent() && !goodRecoil(); }

  G4bool goodNucleus() const;	// Ensure that fragment is energetically okay

protected:
  void fillRecoil();		// Set recoil parameters from CheckBalance
  G4double deltaM() const;	// Mass difference from current parameters

private:
  G4CascadeCheckBalance* balance;	// Used to do kinematics calculations

  G4double excTolerance;	// Minimum excitation energy, rounds to zero

  G4double inputEkin;			// Available initial kinetic energy

  G4int recoilA;			// Nuclear parameters of recoil
  G4int recoilZ;
  G4LorentzVector recoilMomentum;
  G4double excitationEnergy;

  G4ExitonConfiguration theExcitons;	// Used by G4InuclNuclei and G4Fragment

  G4InuclNuclei theRecoilNuclei;	// Reusable buffers for recoil
  G4Fragment theRecoilFragment;

private:
  // Copying of modules is forbidden
  G4CascadeRecoilMaker(const G4CascadeRecoilMaker&);
  G4CascadeRecoilMaker& operator=(const G4CascadeRecoilMaker&);
};

#endif	/* G4CASCADE_RECOIL_MAKER_HH */
