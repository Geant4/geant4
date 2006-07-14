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
// $Id: A01ParaPhysics.cc,v 1.1 2006-07-14 14:43:26 asaim Exp $
// --------------------------------------------------------------
//
// 22-Nov-2004 Construt ALL Particles by T. Koi


#include "A01ParaPhysics.hh"

#include "globals.hh"


A01ParaPhysics::A01ParaPhysics(const G4String& name)
                     :  G4VPhysicsConstructor(name)
{
}

A01ParaPhysics::~A01ParaPhysics()
{
}

void A01ParaPhysics::ConstructParticle()
{
}

#include "G4ParallelWorldScoringProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

void A01ParaPhysics::ConstructProcess()
{
  // Add parallel world scoring process
  G4ParallelWorldScoringProcess* theParallelWorldScoringProcess 
      = new G4ParallelWorldScoringProcess("ParaWorldScoringProc");
  theParallelWorldScoringProcess->SetParallelWorld("ParallelScoringWorld");

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    if (!particle->IsShortLived()) {
      G4ProcessManager* pmanager = particle->GetProcessManager();
      // both postStep and alongStep action are required: because
      // of the use of ghost volumes.
      pmanager->AddProcess(theParallelWorldScoringProcess);
      pmanager->SetProcessOrderingToLast(theParallelWorldScoringProcess, idxAtRest);
      pmanager->SetProcessOrdering(theParallelWorldScoringProcess, idxAlongStep, 1);
      pmanager->SetProcessOrderingToLast(theParallelWorldScoringProcess, idxPostStep);
    }
  }
}


