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
#ifndef G4LORENTZ_CONVERTOR_HH
#define G4LORENTZ_CONVERTOR_HH

#ifndef GLOB
#include "globals.hh"
#endif

#include <vector>

class G4LorentzConvertor {

public:

  G4LorentzConvertor();

  G4LorentzConvertor(const std::vector<G4double>& bmom, 
		     G4double bmass, 
		     const std::vector<G4double>& tmom, 
		     G4double tmass) {

    setBullet(bmom, bmass);
    setTarget(tmom, tmass);
    degenerated = false;  
  }; 

  void setBullet(const std::vector<G4double>& bmom, 
		 G4double bmass) {

    bullet_mom = bmom;
    bullet_mass = bmass;
    //  G4cout << " bullet: e " << bmom[0] << " mass " << bmass << G4endl;
  };

  void setTarget(const std::vector<G4double>& tmom, 
		 G4double tmass) {

    target_mom = tmom;
    target_mass = tmass;
    //  G4cout << " target: e " << tmom[0] << " mass " << tmass << G4endl;
  };

  void toTheCenterOfMass();
 
  void toTheTargetRestFrame(); 

  std::vector<G4double> backToTheLab(const std::vector<G4double>& mom) const; 

  G4double getKinEnergyInTheTRS() const {

    G4double pv = bullet_mom[1] * target_mom[1] + 
      bullet_mom[2] * target_mom[2] + 
      bullet_mom[3] * target_mom[3];  
    G4double ekin_trf = (target_mom[0] * 
			 bullet_mom[0] - pv) / target_mass - bullet_mass;
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
 
  std::vector<G4double> rotate(const std::vector<G4double> mom) const; 

  std::vector<G4double> rotate(const std::vector<G4double> mom1,
			    const std::vector<G4double> mom) const; 

  G4bool reflectionNeeded() const; 

  G4bool trivial() const { 
    return degenerated; 
  }; 

private: 
G4int verboseLevel;
  std::vector<G4double> bullet_mom;
  G4double bullet_mass;

  std::vector<G4double> target_mom;
  G4double target_mass;

  std::vector<G4double> velocity;

  std::vector<G4double> scm_momentum;

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
