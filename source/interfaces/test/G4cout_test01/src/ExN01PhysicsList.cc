// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN01PhysicsList.cc,v 1.1 1999-01-08 16:32:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "ExN01PhysicsList.hh"
#include "G4ParticleTypes.hh"


ExN01PhysicsList::ExN01PhysicsList()
{;}

ExN01PhysicsList::~ExN01PhysicsList()
{;}

void ExN01PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  G4Geantino::GeantinoDefinition();
}

void ExN01PhysicsList::ConstructProcess()
{
  // Define transportation process

  AddTransportation();
}

void ExN01PhysicsList::SetCuts(G4double cut)
{
  // set a cut value to all particles

  SetCutValueForOthers(cut);
}

