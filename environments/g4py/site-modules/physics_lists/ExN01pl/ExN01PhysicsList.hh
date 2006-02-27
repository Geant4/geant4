// $Id: ExN01PhysicsList.hh,v 1.1 2006-02-27 09:50:13 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   ExN01PhysicsList.hh
//
//                                         2005 Q
// ====================================================================
#ifndef EXN01_PHYSICS_LIST_H
#define EXN01_PHYSICS_LIST_H

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class ExN01PhysicsList: public G4VUserPhysicsList {
public:
  ExN01PhysicsList();
  ~ExN01PhysicsList();
  
protected:
  // Construct particle and physics process
  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();
  
};

#endif
