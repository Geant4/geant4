// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4IonConstructor.cc,v 1.2 1999-12-15 14:51:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------
//	GEANT 4 class implementation file 
//

#include "G4IonConstructor.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
// Nuclei
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4GenericIon.hh"

G4IonConstructor::G4IonConstructor()
{
}

G4IonConstructor::~G4IonConstructor()
{
}


void G4IonConstructor::ConstructParticle()
{
  ConstructLightIons();
}

void G4IonConstructor::ConstructLightIons()
{
  //  nuclei
  G4Alpha::AlphaDefinition();
  G4Deuteron::DeuteronDefinition();
  G4Triton::TritonDefinition();
  G4He3::He3Definition();
  //  generic ion
  G4GenericIon::GenericIonDefinition();
}

