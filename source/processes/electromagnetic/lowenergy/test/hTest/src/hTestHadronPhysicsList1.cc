//---------------------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
//---------------------------------------------------------------------------
//
// ClassName:   hTestHadronPhysicsList1
//  
// Description: Implementation file for "LowEnergy" HadronPhysicsList
//
// Authors:    10.04.01 V.Ivanchenko 
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestHadronPhysicsList1.hh"

// Hadron Processes

#include "G4HadronElasticProcess.hh"

#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4LambdaInelasticProcess.hh"
#include "G4SigmaPlusInelasticProcess.hh"
#include "G4SigmaMinusInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

// Low-energy Models

#include "G4LElastic.hh"

#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LELambdaInelastic.hh"
#include "G4LEAntiLambdaInelastic.hh"
#include "G4LESigmaPlusInelastic.hh"
#include "G4LESigmaMinusInelastic.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"

// High-energy Models

#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HEKaonMinusInelastic.hh"
#include "G4HEProtonInelastic.hh"
#include "G4HENeutronInelastic.hh"
#include "G4HELambdaInelastic.hh"
#include "G4HESigmaPlusInelastic.hh"
#include "G4HESigmaMinusInelastic.hh"

// Stopping processes

#ifdef TRIUMF_STOP_PIMINUS
#include "G4PionMinusAbsorptionAtRest.hh"
#else
#include "G4PiMinusAbsorptionAtRest.hh"
#endif
#ifdef TRIUMF_STOP_KMINUS
#include "G4KaonMinusAbsorption.hh"
#else
#include "G4KaonMinusAbsorptionAtRest.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestHadronPhysicsList1::ConstructProcess()
{
   G4cout << "LowEnergy Hadronic PhysicsList is initilized" << G4endl;
   G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess();
   G4LElastic* theElasticModel = new G4LElastic;
   theElasticProcess->RegisterMe(theElasticModel);

   theParticleIterator->reset();
   while ((*theParticleIterator)()) {

      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
       
      if (particleName == "pi+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4PionPlusInelasticProcess* theInelasticProcess =
                 new   G4PionPlusInelasticProcess("inelastic");  
         G4LEPionPlusInelastic* theLEInelasticModel =
                 new G4LEPionPlusInelastic();  
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEPionPlusInelastic* theHEInelasticModel =
                 new G4HEPionPlusInelastic();  
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }

      else if (particleName == "pi-") {   
         pmanager->AddDiscreteProcess(theElasticProcess);   
         G4PionMinusInelasticProcess* theInelasticProcess =
                 new G4PionMinusInelasticProcess("inelastic");
         G4LEPionMinusInelastic* theLEInelasticModel =
                 new G4LEPionMinusInelastic();
         theInelasticProcess->RegisterMe(theLEInelasticModel);

         G4HEPionMinusInelastic* theHEInelasticModel =
                 new G4HEPionMinusInelastic();
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);

#ifdef TRIUMF_STOP_PIMINUS
         pmanager->AddRestProcess(new G4PionMinusAbsorptionAtRest, ordDefault);
#else
         pmanager->AddRestProcess(new G4PiMinusAbsorptionAtRest(
                   G4String("PiMinusAbsorptionAtRest")), ordDefault);
#endif
      }

      else if (particleName == "kaon+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonPlusInelasticProcess* theInelasticProcess =
                 new G4KaonPlusInelasticProcess(G4String("inelastic"));
         G4LEKaonPlusInelastic* theLEInelasticModel =
                 new G4LEKaonPlusInelastic();
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEKaonPlusInelastic* theHEInelasticModel =
                 new G4HEKaonPlusInelastic();
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }

      else if (particleName == "kaon0S") {
         pmanager->AddDiscreteProcess(theElasticProcess);   
         G4KaonZeroSInelasticProcess* theInelasticProcess =
                 new G4KaonZeroSInelasticProcess(G4String("inelastic"));   
         G4LEKaonZeroSInelastic* theLEInelasticModel =
                 new G4LEKaonZeroSInelastic();   
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEKaonZeroInelastic* theHEInelasticModel = 
                 new G4HEKaonZeroInelastic();
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }

      else if (particleName == "kaon0L") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonZeroLInelasticProcess* theInelasticProcess =
                 new G4KaonZeroLInelasticProcess(G4String("inelastic"));
         G4LEKaonZeroLInelastic* theLEInelasticModel =
                 new G4LEKaonZeroLInelastic();
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEKaonZeroInelastic* theHEInelasticModel =
                 new G4HEKaonZeroInelastic();
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }

      else if (particleName == "kaon-") {
         pmanager->AddDiscreteProcess(theElasticProcess);   
         G4KaonMinusInelasticProcess* theInelasticProcess =   
                 new G4KaonMinusInelasticProcess(G4String("inelastic"));
         G4LEKaonMinusInelastic* theLEInelasticModel =
                 new G4LEKaonMinusInelastic();  
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEKaonMinusInelastic* theHEInelasticModel =
                 new G4HEKaonMinusInelastic();
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
#ifdef TRIUMF_STOP_KMINUS
         pmanager->AddRestProcess(new G4KaonMinusAbsorption, ordDefault);
#else
         pmanager->AddRestProcess(new G4KaonMinusAbsorptionAtRest, ordDefault);
#endif
      }

      else if (particleName == "proton") {
         pmanager->AddDiscreteProcess(theElasticProcess);  
         G4ProtonInelasticProcess* theInelasticProcess =
                 new G4ProtonInelasticProcess(G4String("inelastic"));
         G4LEProtonInelastic* theLEInelasticModel = 
                 new  G4LEProtonInelastic();
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEProtonInelastic* theHEInelasticModel = 
                 new G4HEProtonInelastic();
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }

      else if (particleName == "neutron") {
         pmanager->AddDiscreteProcess(theElasticProcess);   
         G4NeutronInelasticProcess* theInelasticProcess =
                 new G4NeutronInelasticProcess(G4String("inelastic"));
         G4LENeutronInelastic* theLEInelasticModel =
                 new G4LENeutronInelastic();
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HENeutronInelastic* theHEInelasticModel =
                 new G4HENeutronInelastic();
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }

      else if (particleName == "lambda") {
         pmanager->AddDiscreteProcess(theElasticProcess);   
         G4LambdaInelasticProcess* theInelasticProcess =
                 new G4LambdaInelasticProcess(G4String("inelastic"));
         G4LELambdaInelastic* theLEInelasticModel = 
                 new G4LELambdaInelastic();
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HELambdaInelastic* theHEInelasticModel = 
                 new  G4HELambdaInelastic();
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }

      else if (particleName == "sigma+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4SigmaPlusInelasticProcess* theInelasticProcess =
                 new G4SigmaPlusInelasticProcess(G4String("inelastic"));
         G4LESigmaPlusInelastic* theLEInelasticModel =
                 new G4LESigmaPlusInelastic();
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HESigmaPlusInelastic* theHEInelasticModel =
                 new G4HESigmaPlusInelastic();  
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }  

      else if (particleName == "sigma-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4SigmaMinusInelasticProcess* theInelasticProcess =
                 new G4SigmaMinusInelasticProcess(G4String("inelastic"));
         G4LESigmaPlusInelastic* theLEInelasticModel =
                 new G4LESigmaPlusInelastic();
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HESigmaMinusInelastic* theHEInelasticModel =
                 new G4HESigmaMinusInelastic();
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);  
      }   
      else if (particleName == "deuteron") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4DeuteronInelasticProcess* theInelasticProcess = 
                            new G4DeuteronInelasticProcess("inelastic");
         G4LEDeuteronInelastic* theLEInelasticModel = 
                                 new G4LEDeuteronInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "triton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4TritonInelasticProcess* theInelasticProcess = 
                            new G4TritonInelasticProcess("inelastic");
         G4LETritonInelastic* theLEInelasticModel = 
                                 new G4LETritonInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "alpha") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AlphaInelasticProcess* theInelasticProcess = 
                            new G4AlphaInelasticProcess("inelastic");
         G4LEAlphaInelastic* theLEInelasticModel = 
                                 new G4LEAlphaInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
    }  
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
