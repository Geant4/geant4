// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst02PhysicsList.cc,v 1.6 2001-01-06 07:04:26 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "Tst02PhysicsList.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include "g4std/iomanip"   

#include "Tst02GeneralPhysics.hh"
#include "Tst02EMPhysics.hh"
#include "Tst02MuonPhysics.hh"
#include "Tst02HadronPhysics.hh"
#include "Tst02IonPhysics.hh"

Tst02PhysicsList::Tst02PhysicsList():  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);

  // General Physics
  RegisterPhysics( new Tst02GeneralPhysics("general") );

  // EM Physics
  RegisterPhysics( new Tst02EMPhysics("standard EM"));

  // Muon Physics
  RegisterPhysics(  new Tst02MuonPhysics("muon"));

   // Hadron Physics
  RegisterPhysics(  new Tst02HadronPhysics("hadron"));

  // Ion Physics
  RegisterPhysics( new Tst02IonPhysics("ion"));

}

Tst02PhysicsList::~Tst02PhysicsList()
{
}

void Tst02PhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}



