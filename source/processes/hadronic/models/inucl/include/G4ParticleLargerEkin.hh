#ifndef G4PARTICLE_LARGER_EKIN_HH
#define G4PARTICLE_LARGER_EKIN_HH

#include "G4InuclElementaryParticle.hh"

class G4ParticleLargerEkin {

public:
  
  G4bool operator() (const G4InuclElementaryParticle& part1,
		     const G4InuclElementaryParticle& part2) {

    return part1.getKineticEnergy() >= part2.getKineticEnergy();
    //  return part1.getEnergy() >= part2.getEnergy();
    //  return part1.getMomModule() >= part2.getMomModule();
  };
 
};

#endif G4PARTICLE_LARGER_EKIN_HH
