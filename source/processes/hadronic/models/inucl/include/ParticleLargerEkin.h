#ifndef PARTICLE_LARGER_EKIN_H
#define PARTICLE_LARGER_EKIN_H

#include "InuclElementaryParticle.h"

class ParticleLargerEkin {

public:
  
bool operator()(const InuclElementaryParticle& part1,
                 const InuclElementaryParticle& part2) {
  return part1.getKineticEnergy() >= part2.getKineticEnergy();
//  return part1.getEnergy() >= part2.getEnergy();
//  return part1.getMomModule() >= part2.getMomModule();
};
 
};

#endif PARTICLE_LARGER_EKIN_H
