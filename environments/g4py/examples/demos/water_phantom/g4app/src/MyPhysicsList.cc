// $Id: MyPhysicsList.cc,v 1.1 2006-02-27 09:44:35 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   MyPhysicsList.cc
//
//                                         2005 Q
// ====================================================================
#include "MyPhysicsList.hh"
#include "Particles.hh"
#include "PhysicsListEMstd.hh"
#include "PhysicsListLHad.hh"

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////////
MyPhysicsList::MyPhysicsList()
  :  G4VModularPhysicsList()
//////////////////////////////
{
  defaultCutValue = 1.*mm;
  SetVerboseLevel(1);

  RegisterPhysics(new Particles);
  RegisterPhysics(new PhysicsListEMstd);
  RegisterPhysics(new PhysicsListLHad);
}

///////////////////////////////
MyPhysicsList::~MyPhysicsList()
///////////////////////////////
{
}

/////////////////////////////
void MyPhysicsList::SetCuts()
/////////////////////////////
{
  SetCutsWithDefault();   
}

