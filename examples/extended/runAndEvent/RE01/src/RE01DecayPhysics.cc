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
//
//
// $Id: RE01DecayPhysics.cc,v 1.1 2004/11/26 07:37:41 asaim Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//

#include "RE01DecayPhysics.hh"

// #include "globals.hh"
// #include "G4ios.hh"
// #include <iomanip>   

RE01DecayPhysics::RE01DecayPhysics(const G4String& name)
                     :  G4VPhysicsConstructor(name)
{
}

RE01DecayPhysics::~RE01DecayPhysics()
{
}

// #include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4UnknownParticle.hh"

void RE01DecayPhysics::ConstructParticle()
{
  G4UnknownParticle::UnknownParticleDefinition();
}

void RE01DecayPhysics::ConstructProcess()
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
    if(particle->GetParticleName()=="unknown") {
      pmanager ->AddProcess(&fUnknownDecay);
      pmanager ->SetProcessOrdering(&fUnknownDecay, idxPostStep);
    }
  }
}


