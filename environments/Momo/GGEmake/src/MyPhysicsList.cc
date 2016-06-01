// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyPhysicsList.cc,v 1.1 1998/09/02 12:10:38 masayasu Exp $
// GEANT4 tag $Name: geant4-00 $
//
// 

#include "MyPhysicsList.hh"
#include "G4ParticleTypes.hh"


MyPhysicsList::MyPhysicsList()
{;}

MyPhysicsList::~MyPhysicsList()
{;}

void MyPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  G4Geantino::GeantinoDefinition();
}

void MyPhysicsList::ConstructProcess()
{
  // Define transportation process

  AddTransportation();
}

void MyPhysicsList::SetCuts(G4double cut)
{
  // set a cut value to all particles

  SetCutValueForOthers(cut);
}


