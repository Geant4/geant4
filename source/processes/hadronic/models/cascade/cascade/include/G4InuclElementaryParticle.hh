#ifndef G4INUCL_ELEMENTARY_PARTICLE_HH
#define G4INUCL_ELEMENTARY_PARTICLE_HH


#include "globals.hh"

#ifndef G4INUCL_PARTICLE_HH
#include "G4InuclParticle.hh"
#endif

class G4InuclElementaryParticle : public G4InuclParticle {

  //      known particle types:
  //      1 - proton
  //      2 - neutron
  //      3 - pi+
  //      5 - pi-
  //      7 - pi 0
  //      10 - photon
  //      111 - quasideutron PP
  //      112 - quasideutron PN
  //      122 - quasideutron NN
 
public:

  G4InuclElementaryParticle() { 

    valid_particle = false;
  };

  G4InuclElementaryParticle(G4int type) 
    : particleType(type) {

    particleMass = getParticleMass(type);
    valid_particle = false;
  };

  G4InuclElementaryParticle(const G4std::vector<G4double>& mom,
			    G4int type) 
    : G4InuclParticle(mom),
      particleType(type) {

    particleMass = getParticleMass(type);
    momentum[0] = sqrt(momentum[1] * momentum[1] + momentum[2] * momentum[2] +
		       momentum[3] * momentum[3] + particleMass * particleMass);
    valid_particle = true;
  };

  G4InuclElementaryParticle(G4double ekin, 
			    G4int type) 
    : particleType(type) {

    particleMass = getParticleMass(type);
    momentum.resize(4);
    momentum[0] = ekin + particleMass;
    momentum[3] = sqrt(momentum[0] * momentum[0] - particleMass * particleMass); 
    momentum[1] = momentum[2] = 0.0;
    valid_particle = true;
  };

  void setType(G4int ityp) { 

    particleType = ityp;
    particleMass = getParticleMass(ityp);
  };

  void setMomentum(const G4std::vector<G4double>& mom) {

    momentum = mom;
    momentum[0] = sqrt(momentum[1] * momentum[1] + momentum[2] * momentum[2] +
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

private: 

  G4int particleType;

  G4double particleMass;

  G4bool valid_particle;

};        

#endif // G4INUCL_ELEMENTARY_PARTICLE_HH 
