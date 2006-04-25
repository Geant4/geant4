// $Id: PhysicsListEMstd.hh,v 1.2 2006-04-25 10:31:40 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   PhysicsListEMstd.hh
//
//   Physics list for electron/positron/gamma
//   EM-standard package
//
//                                         2006 Q
// ====================================================================
#ifndef PHYSICS_LIST_EMSTD_H
#define PHYSICS_LIST_EMSTD_H

#include "G4VUserPhysicsList.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class PhysicsListEMstd: public G4VUserPhysicsList {
public:
  PhysicsListEMstd();
  ~PhysicsListEMstd();

  virtual void ConstructParticle();
  virtual void ConstructProcess();
  virtual void SetCuts();

};

#endif
