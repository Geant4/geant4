// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst01PhysicsList.cc,v 1.1 2001-02-08 08:41:50 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "Tst01PhysicsList.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4ios.hh"
#include "g4std/iomanip"   

#include "Tst01GeneralPhysics.hh"

Tst01PhysicsList::Tst01PhysicsList():  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);

  // General Physics
  RegisterPhysics( new Tst01GeneralPhysics("general") );

}

Tst01PhysicsList::~Tst01PhysicsList()
{
}

void Tst01PhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}



