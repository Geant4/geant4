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
//---------------------------------------------------------------------------
//
//---------------------------------------------------------------------------
//
// ClassName:   hTestHadronPhysicsList1
//  
// Description: Implementation file for "LowEnergy" HadronPhysicsList
//
// Authors:    10.04.01 V.Ivanchenko 
//
// Modified:   11.04.02 V.Ivanchenko migrate to CHIPS
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestHadronPhysicsList1.hh"
#include "G4HadronElasticProcess.hh"
#include "G4LElastic.hh"
#include "G4PionMinusNuclearReaction.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4PionMinusNuclearAtRestChips.hh"
#include "G4KaonMinusAbsorptionAtRest.hh"
#include "G4GammaNuclearReaction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHadronPhysicsList1::ConstructProcess()
{
   G4cout << "CHIPS Hadronic PhysicsList is initilized" << G4endl;
   G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess();
   G4LElastic* theElasticModel = new G4LElastic();
   theElasticProcess->RegisterMe(theElasticModel);

   theParticleIterator->reset();
   while ((*theParticleIterator)()) {

      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
      G4String particleType = particle->GetParticleType();
      G4int bNumber = particle->GetBaryonNumber();
       
      if (particleType == "meson" || bNumber != 0) {
         pmanager->AddDiscreteProcess(theElasticProcess);
         theElasticProcess->RegisterMe(theElasticModel);

         G4String nm = particleName + "Inel";
         G4HadronInelasticProcess* theInelasticProcess = 
                               new G4HadronInelasticProcess(nm, particle);
         G4PionMinusNuclearReaction* chips = new G4PionMinusNuclearReaction();
         theInelasticProcess->RegisterMe(chips);
         pmanager->AddDiscreteProcess(theElasticProcess);
         pmanager->AddDiscreteProcess(theInelasticProcess);

         G4cout << "Process <" << nm << "> is added to " << particleName << G4endl;
      }

      if (particleName == "pi-") {   
         pmanager->AddRestProcess(new G4PionMinusNuclearAtRestChips(), ordDefault);
         G4cout << "Process <PionMinusAnnihilationAtRest> is added to " 
                << particleName << G4endl;
      }

      if (particleName == "kaon-") {
         pmanager->AddRestProcess(new G4KaonMinusAbsorptionAtRest(), ordDefault);
         G4cout << "Process <PionMinusAbsorptionAtRest> is added to " 
                << particleName << G4endl;
      }

      if (particleName == "gamma") {
         G4HadronInelasticProcess* theInelasticProcess = 
                               new G4HadronInelasticProcess("gamNucl", particle);
         theInelasticProcess->RegisterMe(new G4GammaNuclearReaction());
         pmanager->AddDiscreteProcess(theInelasticProcess);
         G4cout << "Process <gamNucl> is added to " 
                << particleName << G4endl;

      }

    }  
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
