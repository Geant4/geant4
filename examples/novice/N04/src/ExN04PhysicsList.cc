// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN04PhysicsList.cc,v 1.9 2001-01-06 06:56:18 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "ExN04PhysicsList.hh"

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

#include "ExN04GeneralPhysics.hh"
#include "ExN04EMPhysics.hh"
#include "ExN04MuonPhysics.hh"
#include "ExN04HadronPhysics.hh"
#include "ExN04IonPhysics.hh"

ExN04PhysicsList::ExN04PhysicsList():  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);

  // General Physics
  RegisterPhysics( new ExN04GeneralPhysics("general") );

  // EM Physics
  RegisterPhysics( new ExN04EMPhysics("standard EM"));

  // Muon Physics
  RegisterPhysics(  new ExN04MuonPhysics("muon"));

   // Hadron Physics
  RegisterPhysics(  new ExN04HadronPhysics("hadron"));

  // Ion Physics
  RegisterPhysics( new ExN04IonPhysics("ion"));


}

ExN04PhysicsList::~ExN04PhysicsList()
{
}

void ExN04PhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}



