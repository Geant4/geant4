#ifndef INUCL_ELEMENTARY_PARTICLE_H
#define INUCL_ELEMENTARY_PARTICLE_H

#include "InuclParticle.h"

class InuclElementaryParticle : public InuclParticle {

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

InuclElementaryParticle() { valid_particle = false;};

InuclElementaryParticle(int type) : particleType(type) {
   particleMass = getParticleMass(type);
   valid_particle = false;
};

InuclElementaryParticle(const vector<double>& mom, int type) :
    InuclParticle(mom), particleType(type) {
   particleMass = getParticleMass(type);
   momentum[0] = sqrt(momentum[1]*momentum[1] + momentum[2]*momentum[2] +
     momentum[3]*momentum[3] + particleMass*particleMass);
   valid_particle = true;
};

InuclElementaryParticle(double ekin, int type) : particleType(type) {
   particleMass = getParticleMass(type);
   momentum.resize(4);
   momentum[0] = ekin + particleMass;
   momentum[3] = sqrt(momentum[0]*momentum[0] - particleMass*particleMass); 
   momentum[1] = momentum[2] = 0.;
   valid_particle = true;
};

void setType(int ityp) { 
  particleType = ityp;
  particleMass = getParticleMass(ityp);
};

void setMomentum(const vector<double>& mom) {
  momentum = mom;
  momentum[0] = sqrt(momentum[1]*momentum[1] + momentum[2]*momentum[2] +
     momentum[3]*momentum[3] + particleMass*particleMass);
  valid_particle = true;
};

int type() const { return particleType; };

bool photon() const { return particleType == 10; };

bool nucleon() const { return particleType <= 2; };

bool pion() const { return particleType == 3 || particleType == 5 
       || particleType == 7; };

bool quasi_deutron() const { return particleType > 100; };

double getMass() const { return particleMass; };

double getParticleMass() const {
  double mass;
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
      mass = 0.;
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
     cout << " uups, unknown particle type " << particleType << endl;
     mass = 0.;
    };
        
    return mass;
};

double getCharge() const {
  double charge;
  switch(particleType) {
    case 1: // proton
      charge = 1.;
      break;
    case 2: // neutron
      charge = 0.;
      break;
    case 3: // pi+
      charge = 1.;
      break;
    case 5: // pi-
      charge = -1.;
      break;
    case 7: // pi0
      charge = 0.;
      break;
    case 10: // photon
      charge = 0.;
      break;
    case 111: // PP
      charge = 2.;
      break;
    case 112: // PN
      charge = 1.;
      break;
    case 122: // NN
      charge = 0.;
      break;
    default:
     cout << " uups, unknown particle type " << particleType << endl;
     charge = 0.;
    };
        
    return charge;
};

double getParticleMass(int type) const {
  double mass;
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
      mass = 0.;
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
     cout << " uups, unknown particle type " << type << endl;
     mass = 0.;
    };
        
    return mass;
};

double getKineticEnergy() const { return momentum[0] - particleMass; };

double getEnergy() const { return momentum[0]; };

bool valid() const { return valid_particle; };

virtual void printParticle() const {
  InuclParticle::printParticle();
  cout << " type " << particleType << " mass " << particleMass << 
    " ekin " << getKineticEnergy() << endl; 
};

private: 

int particleType;

double particleMass;

bool valid_particle;

};        

#endif // INUCL_ELEMENTARY_PARTICLE_H 
