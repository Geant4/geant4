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
// $Id: G4LorentzConvertor.hh 66241 2012-12-13 18:34:42Z gunter $
//
// 20100108  Michael Kelsey -- Use G4LorentzVector internally
// 20100120  M. Kelsey -- BUG FIX:  scm_momentum should be G4ThreeVector
// 20100126  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100519  M. Kelsey -- Add interfaces to pass G4InuclParticles directly
// 20100616  M. Kelsey -- Report bullet and target four-momenta when set
// 20100915  M. Kelsey -- Move constructors to .cc file, add initializers
// 20110602  M. Kelsey -- Drop some unnecessary kinematics intermediates

#ifndef G4LORENTZ_CONVERTOR_HH
#define G4LORENTZ_CONVERTOR_HH

#include "globals.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"

class G4InuclParticle;

class G4LorentzConvertor {
public:
  G4LorentzConvertor();

  G4LorentzConvertor(const G4LorentzVector& bmom, G4double bmass, 
		     const G4LorentzVector& tmom, G4double tmass);

  G4LorentzConvertor(const G4InuclParticle* bullet, 
		     const G4InuclParticle* target);

  void setVerbose(G4int vb=0) { verboseLevel = vb; }

  void setBullet(const G4InuclParticle* bullet);
  void setTarget(const G4InuclParticle* target);

  void setBullet(const G4InuclParticle& bullet) { setBullet(&bullet); }
  void setTarget(const G4InuclParticle& target) { setTarget(&target); }

  // Use correct four-vectors as input
  void setBullet(const G4LorentzVector& bmom) {
    bullet_mom = bmom;
    if (verboseLevel > 3) printBullet();
  }

  void setTarget(const G4LorentzVector& bmom) {
    target_mom = bmom;
    if (verboseLevel > 3) printTarget();
  }

  // These functions "repair" input 4-vectors using specified mass
  void setBullet(const G4LorentzVector& bmom, G4double bmass) {
    bullet_mom.setVectM(bmom.vect(), bmass);
    if (verboseLevel > 3) printBullet();
  }
  
  void setTarget(const G4LorentzVector& tmom, G4double tmass) {
    target_mom.setVectM(tmom.vect(), tmass);
    if (verboseLevel > 3) printTarget();
  }

  // Select reference frame for boosts, rotations, etc.
  void toTheCenterOfMass();
  void toTheTargetRestFrame(); 
  void fillKinematics();	// Common calculations after either of above

  G4LorentzVector backToTheLab(const G4LorentzVector& mom) const;

  // Four-vectors of bullet and target in last chosen reference frame
  const G4LorentzVector& getBullet() const { return bullet_mom; }
  const G4LorentzVector& getTarget() const { return target_mom; }
 
  G4double getKinEnergyInTheTRS() const;
  G4double getTotalSCMEnergy() const { return ecm_tot; }
  G4double getSCMMomentum() const { return scm_momentum.rho(); }
  G4double getTRSMomentum() const;

  G4LorentzVector rotate(const G4LorentzVector& mom) const; 

  G4LorentzVector rotate(const G4LorentzVector& mom1,
			 const G4LorentzVector& mom) const; 

  G4bool reflectionNeeded() const; 

  G4bool trivial() const { return degenerated; }

  // Reporting functions for diagnostics
  void printBullet() const;
  void printTarget() const;

private: 
  static const G4double small;

  G4int verboseLevel;
  G4LorentzVector bullet_mom;
  G4LorentzVector target_mom;

  G4LorentzVector scm_momentum;		// CM momentum relative to target/bullet
  G4ThreeVector   scm_direction;	// Unit vector to reduce repeated calcs

  // Buffer variables for doing ::rotate() calculations
  G4ThreeVector velocity;
  G4double v2;
  G4double ecm_tot;
  G4double valong;
  G4bool degenerated;
};        

#endif // G4LORENTZ_CONVERTOR_HH 
