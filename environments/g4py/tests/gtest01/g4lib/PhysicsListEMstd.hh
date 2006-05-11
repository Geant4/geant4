// $Id: PhysicsListEMstd.hh,v 1.1 2006-05-11 03:00:08 kmura Exp $
// ====================================================================
//   PhysicsListEMstd.hh
//
//                                         2004 Q
// ====================================================================
#ifndef PHYSICS_LIST_EM_STD_H
#define PHYSICS_LIST_EM_STD_H

#include "G4VPhysicsConstructor.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class PhysicsListEMstd : public G4VPhysicsConstructor {

public:
  PhysicsListEMstd();
  ~PhysicsListEMstd();

  virtual void ConstructParticle();
  virtual void ConstructProcess();
};

#endif

