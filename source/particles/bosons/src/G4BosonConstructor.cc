// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BosonConstructor.cc,v 1.1 1999-01-07 16:09:48 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementatBoson file 
//

#include "G4BosonConstructor.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
// Bosons
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"

G4BosonConstructor::G4BosonConstructor()
{
}

G4BosonConstructor::~G4BosonConstructor()
{
}


void G4BosonConstructor::ConstructParticle()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}

