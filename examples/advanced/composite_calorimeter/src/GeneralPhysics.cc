//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
#include "G4ios.hh"
#include "g4std/iomanip"   

GeneralPhysics::GeneralPhysics(const G4String& name)
                     :  G4VPhysicsConstructor(name)
{
}

GeneralPhysics::~GeneralPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
// Bosons
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"

void GeneralPhysics::ConstructParticle()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();  
}

void GeneralPhysics::ConstructProcess()
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


// 2002 by J.P. Wellisch
