// $Id: MyPhysicsList.hh,v 1.1 2006-02-27 09:44:35 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   MyPhysicsList.hh
//
//                                         2005 Q
// ====================================================================
#ifndef MY_PHYSICS_LIST_H
#define MY_PHYSICS_LIST_H

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class MyPhysicsList: public G4VModularPhysicsList {
public:
  MyPhysicsList();
  ~MyPhysicsList();
  
  virtual void SetCuts();
};

#endif
