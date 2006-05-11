// $Id: Particles.cc,v 1.1 2006-05-11 03:00:08 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   Particles.cc
//
//   Physics list for defining particles
//
// ====================================================================
#include "Particles.hh"

#include "G4LeptonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

// ====================================================================
//
// class description
//
// ====================================================================

//////////////////////////////////////
Particles::Particles()
  : G4VPhysicsConstructor("Particles")
//////////////////////////////////////
{
}


///////////////////////
Particles::~Particles()
///////////////////////
{
}


///////////////////////////////////
void Particles::ConstructParticle()
///////////////////////////////////
{
  G4LeptonConstructor::ConstructParticle();
  G4BosonConstructor::ConstructParticle();
  G4MesonConstructor::ConstructParticle();
  G4BaryonConstructor::ConstructParticle();
  G4ShortLivedConstructor::ConstructParticle();
  G4IonConstructor::ConstructParticle(); 
}


//////////////////////////////////
void Particles::ConstructProcess()
//////////////////////////////////
{
}

