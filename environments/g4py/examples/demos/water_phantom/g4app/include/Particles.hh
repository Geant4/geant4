// $Id: Particles.hh,v 1.1 2006-02-27 09:44:35 kmura Exp $
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

