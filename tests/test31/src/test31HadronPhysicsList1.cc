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
//---------------------------------------------------------------------------
//
//---------------------------------------------------------------------------
//
// ClassName:   test31HadronPhysicsList1
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

#include "test31HadronPhysicsList1.hh"
#include "G4HadronElasticProcess.hh"
#include "G4LElastic.hh"
#include "G4PionMinusNuclearReaction.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4PionMinusNuclearAtRestChips.hh"
#include "G4KaonMinusAbsorptionAtRest.hh"
#include "G4GammaNuclearReaction.hh"
#include "G4ProtonAntiProtonAtRestChips.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31HadronPhysicsList1::ConstructProcess()
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

      if (particleName == "anti_proton") {   
         pmanager->AddRestProcess(new G4ProtonAntiProtonAtRestChips(), ordDefault);
         G4cout << "Process <G4ProtonAntiProtonAtRestChips> is added to " 
                << particleName << G4endl;
      }

      if (particleName == "anti_neutron") {   
         pmanager->AddRestProcess(new G4AntiNeutronAnnihilationAtRest(), ordDefault);
         G4cout << "Process <G4AntiNeutronAnnihilationAtRest> is added to " 
                << particleName << G4endl;
      }

      if (particleName == "pi-") {   
         pmanager->AddRestProcess(new G4PionMinusNuclearAtRestChips(), ordDefault);
         G4cout << "Process <G4PionMinusNuclearAtRestChips> is added to " 
                << particleName << G4endl;
      }

      if (particleName == "kaon-") {
         pmanager->AddRestProcess(new G4KaonMinusAbsorptionAtRest(), ordDefault);
         G4cout << "Process <G4KaonMinusAbsorptionAtRest> is added to " 
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
