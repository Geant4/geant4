// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst02GeneralPhysics.cc,v 1.1 2001-01-06 07:04:25 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "Tst02GeneralPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/iomanip"   

Tst02GeneralPhysics::Tst02GeneralPhysics(const G4String& name)
                     :  G4VPhysicsConstructor(name)
{
}

Tst02GeneralPhysics::~Tst02GeneralPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
// Bosons
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"

void Tst02GeneralPhysics::ConstructParticle()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();  
}

void Tst02GeneralPhysics::ConstructProcess()
{
  // Add Decay Process
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (fDecayProcess.IsApplicable(*particle)) { 
      pmanager ->AddProcess(&fDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(&fDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(&fDecayProcess, idxAtRest);
    }
  }
}


