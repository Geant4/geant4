// $Id: QPhysicsList.cc,v 1.1 2006-02-27 09:50:13 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   QPhysicsList.cc
//
//                                         2005 Q
// ====================================================================
#include "QPhysicsList.hh"
#include "PhysicsListEMstd.hh"

// ====================================================================
//
// class description
//
// ====================================================================

////////////////////////////
QPhysicsList::QPhysicsList()
  :  G4VModularPhysicsList()
////////////////////////////
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.*mm;
  SetVerboseLevel(1);

  // EM Physics
  RegisterPhysics(new PhysicsListEMstd);
}

/////////////////////////////
QPhysicsList::~QPhysicsList()
/////////////////////////////
{
}

////////////////////////////
void QPhysicsList::SetCuts()
////////////////////////////
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}

