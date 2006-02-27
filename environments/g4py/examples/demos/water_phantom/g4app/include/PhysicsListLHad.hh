// $Id: PhysicsListLHad.hh,v 1.1 2006-02-27 09:44:35 kmura Exp $
// ====================================================================
//   PhysicsListLHad.hh
//
//                                         2004 Q
// ====================================================================
#ifndef PHYSICS_L_HAD_H
#define PHYSICS_L_HAD_H

#include "G4VPhysicsConstructor.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class PhysicsListLHad : public G4VPhysicsConstructor {

public:
  PhysicsListLHad();
  ~PhysicsListLHad();

  virtual void ConstructParticle();
  virtual void ConstructProcess();
};

#endif

