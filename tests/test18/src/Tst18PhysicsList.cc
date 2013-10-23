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
// $Id$
//

#include "Tst18PhysicsList.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"

#include "G4IonTable.hh"
#include "G4Ions.hh"

#include "G4RadioactiveDecay.hh"


//#include "G4FastSimulationManagerProcess.hh"


Tst18PhysicsList::Tst18PhysicsList():  G4VUserPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.*m;

  SetVerboseLevel(1);
}

Tst18PhysicsList::~Tst18PhysicsList()
{
}

void Tst18PhysicsList::ConstructParticle()
{

  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  // create all particles
  ConstructAllBosons();
  ConstructAllLeptons();
  ConstructAllMesons();
  ConstructAllBaryons();
  ConstructAllIons();
  ConstructAllShortLiveds();

  G4GenericIon::GenericIonDefinition();
}

void Tst18PhysicsList::ConstructProcess()
{
  AddTransportation();

  ConstructEM();
  ConstructHad();
  ConstructGeneral();
}

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

void Tst18PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma
      // Construct processes for gamma
      //      pmanager->AddDiscreteProcess(new G4GammaConversion());
      //pmanager->AddDiscreteProcess(new G4ComptonScattering());      
      //pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } else if (particleName == "e-") {
    //electron
      pmanager->AddProcess(new G4eMultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,3);

    } else if (particleName == "e+") {
    //positron
      pmanager->AddProcess(new G4eMultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
    //muon  
      pmanager->AddProcess(new G4MuMultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4MuIonisation,      -1, 2,2);
      pmanager->AddProcess(new G4MuBremsstrahlung,  -1,-1,3);
      pmanager->AddProcess(new G4MuPairProduction,  -1,-1,4);
     
   } else if( particleName == "GenericIon" ) {
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4ionIonisation,      -1, 2,2);

 
   } else if ((!particle->IsShortLived()) &&
	      (particle->GetPDGCharge() != 0.0) && 
	      (particle->GetParticleName() != "chargedgeantino")) {
     // all others charged particles except geantino

      pmanager->AddProcess(new G4hMultipleScattering(),-1,1,1);
      G4hIonisation*    hIon;
      hIon= new G4hIonisation();
      pmanager->AddProcess(hIon,  -1,2,2);      

      G4double demax = 0.05;  // try to lose at most 5% of the energy in 
                              //    a single step (in limit of large energies)
      G4double stmin = 0.01 * mm;  // length of the final step 
      hIon->SetStepFunction( demax, stmin );
    }
  }
}


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
#include "G4XiZeroInelasticProcess.hh"
#include "G4XiMinusInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4OmegaMinusInelasticProcess.hh"

// Models
#include "G4HadronElastic.hh"
#include "G4CascadeInterface.hh"
#include "G4BinaryLightIonReaction.hh"

// Stopping processes
#include "G4PiMinusAbsorptionBertini.hh"
#include "G4KaonMinusAbsorptionBertini.hh"


void Tst18PhysicsList::ConstructHad()
{
  // Elastic model and process
  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  G4HadronElastic* theElasticModel = new G4HadronElastic;
  theElasticProcess->RegisterMe(theElasticModel);

  // Inelastic hadronic model for hadrons
  G4CascadeInterface* bertini = new G4CascadeInterface;

  // Inelastic hadronic model for ions
  G4BinaryLightIonReaction* binaryCascade = new G4BinaryLightIonReaction;
  binaryCascade->SetMinEnergy(0.0);
  binaryCascade->SetMaxEnergy(110*MeV);

  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "pi+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4PionPlusInelasticProcess* theInelasticProcess = 
                             new G4PionPlusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "pi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4PionMinusInelasticProcess* theInelasticProcess = 
                             new G4PionMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4PiMinusAbsorptionBertini, ordDefault);

    } else if (particleName == "kaon+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonPlusInelasticProcess* theInelasticProcess = 
                             new G4KaonPlusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon0S") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonZeroSInelasticProcess* theInelasticProcess = 
                             new G4KaonZeroSInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon0L") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonZeroLInelasticProcess* theInelasticProcess = 
                             new G4KaonZeroLInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonMinusInelasticProcess* theInelasticProcess = 
                             new G4KaonMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4KaonMinusAbsorptionBertini, ordDefault);

    } else if (particleName == "proton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4ProtonInelasticProcess* theInelasticProcess = 
                             new G4ProtonInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "neutron") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4NeutronInelasticProcess* theInelasticProcess = 
                             new G4NeutronInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "lambda") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LambdaInelasticProcess* theInelasticProcess = 
                             new G4LambdaInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "sigma+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4SigmaPlusInelasticProcess* theInelasticProcess = 
                             new G4SigmaPlusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "sigma-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4SigmaMinusInelasticProcess* theInelasticProcess = 
                             new G4SigmaMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "xi0") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4XiZeroInelasticProcess* theInelasticProcess = 
                             new G4XiZeroInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "xi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4XiMinusInelasticProcess* theInelasticProcess = 
                             new G4XiMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "omega-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4OmegaMinusInelasticProcess* theInelasticProcess =
                            new G4OmegaMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "deuteron") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4DeuteronInelasticProcess* theInelasticProcess = 
                            new G4DeuteronInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(binaryCascade);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "triton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4TritonInelasticProcess* theInelasticProcess = 
                            new G4TritonInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(binaryCascade);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "alpha") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AlphaInelasticProcess* theInelasticProcess = 
                            new G4AlphaInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(binaryCascade);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
  }
}


#include "G4Decay.hh"
void Tst18PhysicsList::ConstructGeneral()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
  // Declare radioactive decay to the GenericIon in the IonTable.
  //
/////////  const G4IonTable *theIonTable = G4ParticleTable::GetParticleTable()->GetIonTable();
/////////  for (G4int i = 0; i < theIonTable->Entries(); i++) {
/////////    G4String particleName = theIonTable->GetParticle(i)->GetParticleName();
/////////    if (particleName == "GenericIon") {
  G4ParticleDefinition* gion = G4ParticleTable::GetParticleTable()->GetGenericIon();
  if(gion)
  {
      G4RadioactiveDecay *theRadioactiveDecay = new G4RadioactiveDecay();
      G4ProcessManager* pmanager = gion->GetProcessManager();
      pmanager->SetVerboseLevel(0);
      pmanager ->AddProcess(theRadioactiveDecay);
      pmanager ->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
      pmanager ->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
  }
/////////    }
/////////  } 
}

void Tst18PhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}

void Tst18PhysicsList::ConstructAllBosons()
{
  // Construct all bosons
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst18PhysicsList::ConstructAllLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst18PhysicsList::ConstructAllMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst18PhysicsList::ConstructAllBaryons()
{
  //  Construct all barions
  G4BaryonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst18PhysicsList::ConstructAllIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void Tst18PhysicsList::ConstructAllShortLiveds()
{
  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

