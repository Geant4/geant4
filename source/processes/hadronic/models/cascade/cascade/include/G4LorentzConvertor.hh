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
// $Id: G4LorentzConvertor.hh,v 1.11 2010-01-22 00:07:10 mkelsey Exp $
//
// 20100108  Michael Kelsey -- Use G4LorentzVector internally
// 20100120  M. Kelsey -- BUG FIX:  scm_momentum should be G4ThreeVector

#ifndef G4LORENTZ_CONVERTOR_HH
#define G4LORENTZ_CONVERTOR_HH

#ifndef GLOB
#include "globals.hh"
#endif

#include <vector>
#include "G4CascadeMomentum.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"


class G4LorentzConvertor {

public:

  G4LorentzConvertor();

  G4LorentzConvertor(const G4CascadeMomentum& bmom, 
		     G4double bmass, 
		     const G4CascadeMomentum& tmom, 
		     G4double tmass) {
    setBullet(bmom, bmass);
    setTarget(tmom, tmass);
    degenerated = false;  
  }; 

  void setBullet(const G4CascadeMomentum& bmom, 
		 G4double bmass) {
    bullet_mom.setVectM(G4ThreeVector(bmom[1],bmom[2],bmom[3]), bmass);
    //  G4cout << " bullet: e " << bullet_mom.e() << " mass "
    //         << bullet_mom.m() << G4endl;
  };

  void setTarget(const G4CascadeMomentum& tmom, 
		 G4double tmass) {
    target_mom.setVectM(G4ThreeVector(tmom[1],tmom[2],tmom[3]), tmass);
    //  G4cout << " target: e " << target_mom.e() << " mass "
    //         << target_mom.m() << G4endl;
  };

  void toTheCenterOfMass();
 
  void toTheTargetRestFrame(); 

  G4CascadeMomentum backToTheLab(const G4LorentzVector& mom) const; 

  G4double getKinEnergyInTheTRS() const {
    G4double pv = bullet_mom.vect().dot(target_mom.vect());
    G4double ekin_trf = (target_mom.e() * bullet_mom.e() - pv)
                          / target_mom.m() - bullet_mom.m();
    return ekin_trf; 
  };

  G4double getTotalSCMEnergy() const { 
    return ecm_tot; 
  };

  G4double getSCMMomentum() const { 
    return pscm; 
  };

  G4double getTRSMomentum() const { 
    return plab; 
  };
 
  G4CascadeMomentum rotate(const G4LorentzVector& mom) const; 

  G4CascadeMomentum rotate(const G4LorentzVector& mom1,
			    const G4LorentzVector& mom) const; 

  G4bool reflectionNeeded() const; 

  G4bool trivial() const { 
    return degenerated; 
  }; 

private: 
  static const G4double small;

  G4int verboseLevel;
  G4LorentzVector bullet_mom;
  G4LorentzVector target_mom;

  G4ThreeVector velocity;

  G4ThreeVector scm_momentum;

  G4double ecm_tot;

  G4double pscm;

  G4double plab;

  G4double gamma;

  G4double v2;

  G4double ga;

  G4double gb;

  G4double gbpp;

  G4double gapp;

  G4bool degenerated;
};        

#endif // G4LORENTZ_CONVERTOR_HH 
