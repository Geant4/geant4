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
/// \file field/field06/src/F06ExtraPhysics.cc
/// \brief Implementation of the F06ExtraPhysics class
//
//
//

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"

#include "G4Decay.hh"
#include "G4DecayTable.hh"
#include "G4NeutronBetaDecayChannel.hh"

#include "G4UserSpecialCuts.hh"
#include "G4StepLimiter.hh"
#include "G4SystemOfUnits.hh"

#include "F06ExtraPhysics.hh"

F06ExtraPhysics::F06ExtraPhysics() 
    : G4VPhysicsConstructor("Extra") {;}

F06ExtraPhysics::~F06ExtraPhysics() {;}

void F06ExtraPhysics::ConstructParticle() {;}

void F06ExtraPhysics::ConstructProcess()
{
    theParticleIterator->reset();

    while ((*theParticleIterator)()) {
        G4ParticleDefinition* particle = theParticleIterator->value();
        G4ProcessManager* pmanager = particle->GetProcessManager();
        G4String particleName = particle->GetParticleName();

        if (!pmanager) {
            std::ostringstream o;
            o << "Particle " << particleName << "without a Process Manager";
            G4Exception("F06ExtraPhysics::ConstructProcess()","",
                         FatalException,o.str().c_str());
        }

        pmanager->AddDiscreteProcess(new G4StepLimiter());
        pmanager->AddDiscreteProcess(new G4UserSpecialCuts());
    }

    AddBetaDecay();
}

void F06ExtraPhysics::AddBetaDecay()
{
    theParticleIterator->reset();

    while ((*theParticleIterator)()) {

        G4ParticleDefinition* particle = theParticleIterator->value();
        G4String particleName = particle->GetParticleName();

        if (particleName == "neutron") {

           particle->SetPDGLifeTime(885.7*second);
           particle->SetPDGStable(false);

           G4DecayTable * table = new G4DecayTable();
           G4VDecayChannel* mode = new G4NeutronBetaDecayChannel("neutron",1.00);
           table->Insert(mode);
           particle->SetDecayTable(table);

           G4ProcessManager* pmanager = particle->GetProcessManager();
           if (!pmanager) {
               std::ostringstream o;
               o << "Particle " << particleName << "without a Process Manager";
               G4Exception("F06ExtraPhysics::ConstructProcess()","",
                            FatalException,o.str().c_str());
           }

           G4Decay* theDecayProcess = new G4Decay();
           pmanager->AddProcess(theDecayProcess);
           pmanager->SetProcessOrdering(theDecayProcess,idxPostStep);
           pmanager->SetProcessOrdering(theDecayProcess,idxAtRest);

           break;
        }
    }
}
