// $Id: Particles.hh,v 1.1 2006-05-11 04:35:32 kmura Exp $
// ====================================================================
//   Particles.hh
//
//                                         2004 Q
// ====================================================================
#ifndef PARTICLES_H
#define PARTICLES_H

#include "G4VPhysicsConstructor.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class Particles : public G4VPhysicsConstructor {

public:
  Particles();
  ~Particles();

  virtual void ConstructParticle();
  virtual void ConstructProcess();
};

#endif

