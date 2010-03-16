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
// $Id: G4LorentzConvertor.hh,v 1.14 2010-03-16 22:10:26 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100108  Michael Kelsey -- Use G4LorentzVector internally
// 20100120  M. Kelsey -- BUG FIX:  scm_momentum should be G4ThreeVector
// 20100126  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly

#ifndef G4LORENTZ_CONVERTOR_HH
#define G4LORENTZ_CONVERTOR_HH

#ifndef GLOB
#include "globals.hh"
#endif

#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"

class G4LorentzConvertor {

public:

  G4LorentzConvertor();

  G4LorentzConvertor(const G4LorentzVector& bmom, G4double bmass, 
		     const G4LorentzVector& tmom, G4double tmass) {
    setBullet(bmom, bmass);
    setTarget(tmom, tmass);
  }; 

  void setVerbose(G4int vb=0) { verboseLevel = vb; }

  // NOTE:  These functions "repair" input 4-vectors using specified mass
  void setBullet(const G4LorentzVector& bmom, G4double bmass) {
    bullet_mom.setVectM(bmom.vect(), bmass);

    //  G4cout << " bullet: e " << bullet_mom.e() << " mass "
    //         << bullet_mom.m() << G4endl;
  };

  void setTarget(const G4LorentzVector& tmom, G4double tmass) {
    target_mom.setVectM(tmom.vect(), tmass);

    //  G4cout << " target: e " << target_mom.e() << " mass "
    //         << target_mom.m() << G4endl;
  };

  void toTheCenterOfMass();
  void toTheTargetRestFrame(); 

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

private: 
  static const G4double small;

  G4int verboseLevel;
  G4LorentzVector bullet_mom;
  G4LorentzVector target_mom;

  G4LorentzVector scm_momentum;		// CM momentum relative to target/bullet

  // Buffer variables for doing ::rotate() calculations
  G4ThreeVector velocity;
  G4double gamma;
  G4double v2;
  G4double ecm_tot;
  G4double ga;
  G4double gb;
  G4double gbpp;
  G4double gapp;
  G4bool degenerated;
};        

#endif // G4LORENTZ_CONVERTOR_HH 
