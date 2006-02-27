// $Id: QPhysicsList.hh,v 1.1 2006-02-27 09:50:13 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   QPhysicsList.hh
//
//                                         2005 Q
// ====================================================================
#ifndef Q_PHYSICS_LIST_H
#define Q_PHYSICS_LIST_H

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

// ====================================================================
//
// class definition
//
// ====================================================================

class QPhysicsList: public G4VModularPhysicsList {
public:
  QPhysicsList();
  ~QPhysicsList();
  
  virtual void SetCuts();
};

#endif
