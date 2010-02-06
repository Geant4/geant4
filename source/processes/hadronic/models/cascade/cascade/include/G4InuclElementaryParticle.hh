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
#ifndef G4INUCL_ELEMENTARY_PARTICLE_HH
#define G4INUCL_ELEMENTARY_PARTICLE_HH


#include "globals.hh"

#ifndef G4INUCL_PARTICLE_HH
#include "G4InuclParticle.hh"
#endif

class G4InuclElementaryParticle : public G4InuclParticle {

//                     known particle types:
//      1 - proton          11 - k+         111 - quasideuteron PP
//      2 - neutron         13 - k-         112 - quasideuteron PN
//      3 - pi+             15 - k0         122 - quasideuteron NN
//      5 - pi-             17 - k0bar
//      7 - pi 0            21 - lambda 
//     10 - photon          23 - sigma+
//                          25 - sigma0
//                          27 - sigma-
//                          29 - xi0
//                          31 - xi-
 
public:

  G4InuclElementaryParticle() { 
    particleType = 0;     // DHW: added to keep 4.3 compiler happy
    particleMass = 0.;    //            "              "
    valid_particle = false;
    generation = 0;
  };

  G4InuclElementaryParticle(G4int type) 
    : particleType(type) {

    particleMass = getParticleMass(type);
    valid_particle = false;
  };

  G4InuclElementaryParticle(const G4CascadeMomentum& mom,
			    G4int type) 
    : G4InuclParticle(mom),
      particleType(type) {

    particleMass = getParticleMass(type);
    momentum[0] = std::sqrt(momentum[1] * momentum[1] + momentum[2] * momentum[2] +
		       momentum[3] * momentum[3] + particleMass * particleMass);
    valid_particle = true;
  };

  
  G4InuclElementaryParticle(const G4CascadeMomentum& mom,
			    G4int type, G4int model) 
    : G4InuclParticle(mom),
      particleType(type) {

    G4InuclParticle::setModel(model);

    particleMass = getParticleMass(type);
    momentum[0] = std::sqrt(momentum[1] * momentum[1] + momentum[2] * momentum[2] +
		       momentum[3] * momentum[3] + particleMass * particleMass);
    valid_particle = true;
  };

  G4InuclElementaryParticle(G4double ekin, 
			    G4int type) 
    : particleType(type) {

    particleMass = getParticleMass(type);
    momentum[0] = ekin + particleMass;
    momentum[3] = std::sqrt(momentum[0] * momentum[0] - particleMass * particleMass); 
    momentum[1] = momentum[2] = 0.0;
    valid_particle = true;
  };

  void setType(G4int ityp) { 

    particleType = ityp;
    particleMass = getParticleMass(ityp);
  };

  void setMomentum(const G4CascadeMomentum& mom) {

    momentum = mom;
    momentum[0] = std::sqrt(momentum[1] * momentum[1] + momentum[2] * momentum[2] +
		       momentum[3] * momentum[3] + particleMass * particleMass);
    valid_particle = true;
  };

  G4int type() const { 

    return particleType; 
  };

  G4bool photon() const { 

    return particleType == 10; 
  };

  G4bool nucleon() const { 
    return particleType <= 2; 
  };

  G4bool baryon() const { 
    return (particleType == 1  ||
            particleType == 2  ||
            particleType == 21 ||
            particleType == 23 ||
            particleType == 25 ||
            particleType == 27 ||
            particleType == 29 ||
            particleType == 31 );
  };

  G4bool pion() const { 

    return particleType == 3 || particleType == 5 || particleType == 7; 
  };

  G4bool quasi_deutron() const { 

    return particleType > 100; 
  };

  G4double getMass() const { 

    return particleMass; 
  };

  G4double getParticleMass() const {

    G4double mass;

    switch(particleType) {
    case 1: // proton
      mass = 0.93827;
      break;
    case 2: // neutron
      mass = 0.93957;
      break;
    case 3: // pi+
      mass = 0.13957;
      break;
    case 5: // pi-
      mass = 0.13957;
      break;
    case 7: // pi0
      mass = 0.13498;
      break;
    case 10: // photon
      mass = 0.0;
      break;
    case 11: // k+
      mass = 0.49368;
      break;
    case 13: // k-
      mass = 0.49368;
      break;
    case 15: // k0
      mass = 0.49767;
      break;
    case 17: // k0bar
      mass = 0.49767;
      break;
    case 21: // lambda
      mass = 1.1157;
      break;
    case 23: // sigma+
      mass = 1.1894;
      break;
    case 25: // sigma0
      mass = 1.1926;
      break;
    case 27: // sigma-
      mass = 1.1974;
      break;
    case 29: // xi0
      mass = 1.3148;
      break;
    case 31: // xi-
      mass = 1.3213;
      break;
    case 111: // PP
      mass = 0.93827 + 0.93827;
      break;
    case 112: // PN
      mass = 0.93827 + 0.93957;
      break;
    case 122: // NN
      mass = 0.93957 + 0.93957;
      break;
    default:
      G4cout << " uups, unknown particle type " << particleType << G4endl;
      mass = 0.;
    };
        
    return mass;
  };

  G4double getCharge() const {

    G4double charge;

    switch(particleType) {
    case 1: // proton
      charge = 1.0;
      break;
    case 2: // neutron
      charge = 0.0;
      break;
    case 3: // pi+
      charge = 1.0;
      break;
    case 5: // pi-
      charge = -1.0;
      break;
    case 7: // pi0
      charge = 0.0;
      break;
    case 10: // photon
      charge = 0.0;
      break;
    case 11: // k+
      charge = 1.0;
      break;
    case 13: // k-
      charge = -1.0;
      break;
    case 15: // k0
      charge = 0.0;
      break;
    case 17: // k0bar
      charge = 0.0;
      break;
    case 21: // lambda
      charge = 0.0;
      break;
    case 23: // sigma+
      charge = 1.0;
      break;
    case 25: // sigma0
      charge = 0.0;
      break;
    case 27: // sigma-
      charge = -1.0;
      break;
    case 29: // xi0
      charge = 0.0;
      break;
    case 31: // xi-
      charge = -1.0;
      break;
    case 111: // PP
      charge = 2.0;
      break;
    case 112: // PN
      charge = 1.0;
      break;
    case 122: // NN
      charge = 0.0;
      break;
    default:
      G4cout << " uups, unknown particle type " << particleType << G4endl;
      charge = 0.0;
    };
        
    return charge;
  };


  G4double getStrangeness(G4int type) const {

    G4double strangeness;

    switch(type) {
    case 1: // proton
      strangeness = 0.0;
      break;
    case 2: // neutron
      strangeness = 0.0;
      break;
    case 3: // pi+
      strangeness = 0.0;
      break;
    case 5: // pi-
      strangeness = 0.0;
      break;
    case 7: // pi0
      strangeness = 0.0;
      break;
    case 10: // photon
      strangeness = 0.0;
      break;
    case 11: // k+
      strangeness = 1.0;
      break;
    case 13: // k-
      strangeness = -1.0;
      break;
    case 15: // k0
      strangeness = 1.0;
      break;
    case 17: // k0bar
      strangeness = -1.0;
      break;
    case 21: // lambda
      strangeness = -1.0;
      break;
    case 23: // sigma+
      strangeness = -1.0;
      break;
    case 25: // sigma0
      strangeness = -1.0;
      break;
    case 27: // sigma-
      strangeness = -1.0;
      break;
    case 29: // xi0
      strangeness = -2.0;
      break;
    case 31: // xi-
      strangeness = -2.0;
      break;
    case 111: // PP
      strangeness = 0.0;
      break;
    case 112: // PN
      strangeness = 0.0;
      break;
    case 122: // NN
      strangeness = 0.0;
      break;
    default:
      G4cout << " unknown particle type " << type << G4endl;
      strangeness = 0.0;
    };
        
    return strangeness;
  };


  G4double getParticleMass(G4int type) const {

    G4double mass;

    switch(type) {
    case 1: // proton
      mass = 0.93827;
      break;
    case 2: // neutron
      mass = 0.93957;
      break;
    case 3: // pi+
      mass = 0.13957;
      break;
    case 5: // pi-
      mass = 0.13957;
      break;
    case 7: // pi0
      mass = 0.13498;
      break;
    case 10: // photon
      mass = 0.0;
      break;
    case 11: // k+
      mass = 0.49368;
      break;
    case 13: // k-
      mass = 0.49368;
      break;
    case 15: // k0
      mass = 0.49767;
      break;
    case 17: // k0bar
      mass = 0.49767;
      break;
    case 21: // lambda
      mass = 1.1157;
      break;
    case 23: // sigma+
      mass = 1.1894;
      break;
    case 25: // sigma0
      mass = 1.1926;
      break;
    case 27: // sigma-
      mass = 1.1974;
      break;
    case 29: // xi0
      mass = 1.3148;
      break;
    case 31: // xi-
      mass = 1.3213;
      break;
    case 111: // PP
      mass = 0.93827 + 0.93827;
      break;
    case 112: // PN
      mass = 0.93827 + 0.93957;
      break;
    case 122: // NN
      mass = 0.93957 + 0.93957;
      break;
    default:
      G4cout << " uups, unknown particle type " << type << G4endl;
      mass = 0.0;
    };
        
    return mass;
  };

  G4double getKineticEnergy() const { 

    return momentum[0] - particleMass; 
  };

  G4double getEnergy() const { 

    return momentum[0]; 
  };

  G4bool valid() const { 

    return valid_particle; 
  };

  virtual void printParticle() const {

    G4InuclParticle::printParticle();

    G4cout << " Particle: type " << particleType << " mass " << particleMass << 
      " ekin " << getKineticEnergy() << G4endl; 
  };

  void setGeneration(G4int gen) {
    generation = gen;
  }

  G4int getGeneration() {
    return generation;
  }

private: 

  G4int particleType;

  G4double particleMass;

  G4bool valid_particle;

  G4int generation;

};        

#endif // G4INUCL_ELEMENTARY_PARTICLE_HH 
