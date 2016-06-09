//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: RE01DecayPhysics.cc,v 1.2 2006-06-29 17:43:43 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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


