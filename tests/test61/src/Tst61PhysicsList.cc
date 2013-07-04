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
// $Id$
//
// 091118 Change multiple scattering processes to particle dedicated by T. Koi
//
#include <iomanip>

#include "globals.hh"
#include "Tst61PhysicsList.hh"
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
#include "G4IonsShenCrossSection.hh"


Tst61PhysicsList::Tst61PhysicsList():  G4VUserPhysicsList()
{
  SetVerboseLevel(1);
}

Tst61PhysicsList::~Tst61PhysicsList()
{
}

void Tst61PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program.

  ConstructAllBosons();
  ConstructAllLeptons();
  ConstructAllMesons();
  ConstructAllBaryons();
  ConstructAllIons();
  ConstructAllShortLiveds();
}

void Tst61PhysicsList::ConstructAllBosons()
{
  // Construct all bosons
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst61PhysicsList::ConstructAllLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst61PhysicsList::ConstructAllMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst61PhysicsList::ConstructAllBaryons()
{
  //  Construct all barions
  G4BaryonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst61PhysicsList::ConstructAllIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst61PhysicsList::ConstructAllShortLiveds()
{
  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst61PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructLeptHad();
  ConstructHad();
  ConstructGeneral();
}

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

//#include "G4MultipleScattering.hh"
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

void Tst61PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      // Construct processes for gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } else if (particleName == "e-") {
      // Construct processes for electron
      pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
      pmanager->AddProcess(new G4eIonisation(),-1,2,2);
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);

    } else if (particleName == "e+") {
      // Construct processes for positron
      pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);

      pmanager->AddProcess(new G4eIonisation(),-1,2,2);
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation(),0,-1,4);

    } else if( particleName == "mu+" ||
               particleName == "mu-"    ) {
     // Construct processes for muon+
      pmanager->AddProcess(new G4MuMultipleScattering(),-1,1,1);
      pmanager->AddProcess(new G4MuIonisation(),-1,2,2);
      pmanager->AddProcess(new G4MuBremsstrahlung(),-1,-1,3);
      pmanager->AddProcess(new G4MuPairProduction(),-1,-1,4);

    } else if( particleName == "GenericIon" ) {
      pmanager->AddProcess(new G4hMultipleScattering(),-1,1,1);
      pmanager->AddProcess(new G4hIonisation(),-1,2,2);

    } else {
      if ((particle->GetPDGCharge() != 0.0) &&
          (particle->GetParticleName() != "chargedgeantino")&&
          (!particle->IsShortLived()) ) {
         // all others charged particles except geantino
        pmanager->AddProcess(new G4hMultipleScattering(),-1,1,1);
        pmanager->AddProcess(new G4hIonisation(),-1,2,2);
      }
    }
  }
}

// Hadron Processes

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4LambdaInelasticProcess.hh"
#include "G4AntiLambdaInelasticProcess.hh"
#include "G4SigmaPlusInelasticProcess.hh"
#include "G4SigmaMinusInelasticProcess.hh"
#include "G4AntiSigmaPlusInelasticProcess.hh"
#include "G4AntiSigmaMinusInelasticProcess.hh"
#include "G4XiZeroInelasticProcess.hh"
#include "G4XiMinusInelasticProcess.hh"
#include "G4AntiXiZeroInelasticProcess.hh"
#include "G4AntiXiMinusInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4IonInelasticProcess.hh"
#include "G4OmegaMinusInelasticProcess.hh"
#include "G4AntiOmegaMinusInelasticProcess.hh"

// Low-energy Models

#include "G4HadronElastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"

#include "G4QMDReaction.hh"

// Generator models
#include "G4Evaporation.hh"
#include "G4CompetitiveFission.hh"
#include "G4FermiBreakUp.hh"
#include "G4StatMF.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4LEProtonInelastic.hh"
#include "G4StringModel.hh"
#include "G4PreCompoundModel.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

// -- bc
//#include "G4BinaryCascade.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4TripathiCrossSection.hh"

#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4CascadeInterface.hh"
#include "G4BinaryLightIonReaction.hh"


void Tst61PhysicsList::ConstructHad()
{
  //Most hadrons
  //Bertini at low energies, then FTFP
  //theFTFPwtBERT + theBertini

  G4TheoFSGenerator* theFTFPwtBERT = new G4TheoFSGenerator("FTFP");
  G4PreCompoundModel* thePreEquilib;
  G4ExcitationHandler* theHandler;
  G4GeneratorPrecompoundInterface* theCascade;
  G4FTFModel* theStringModel;
  G4ExcitedStringDecay* theStringDecay;
  G4LundStringFragmentation* theLund;

  theStringModel = new G4FTFModel;
  theStringDecay = new G4ExcitedStringDecay(theLund = new G4LundStringFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);
 
  theCascade = new G4GeneratorPrecompoundInterface;
  thePreEquilib = new G4PreCompoundModel(theHandler = new G4ExcitationHandler);
  theCascade->SetDeExcitation(thePreEquilib);  
 
  theFTFPwtBERT->SetTransport(theCascade);
  theFTFPwtBERT->SetHighEnergyGenerator(theStringModel);
  theFTFPwtBERT->SetMinEnergy(2.*GeV);
  theFTFPwtBERT->SetMaxEnergy(100.*TeV);
 
  G4CascadeInterface* theBertini = new G4CascadeInterface;
  theBertini->SetMinEnergy(0.*GeV);
  theBertini->SetMaxEnergy(6.*GeV);
 
  //AntiHyperons:
  //Use FTFP for full energy range
  G4TheoFSGenerator* theFTFPanti = new G4TheoFSGenerator("FTFP");
  theFTFPanti->SetMinEnergy(0.*GeV);
  theFTFPanti->SetMaxEnergy(100.*TeV);
  theFTFPanti->SetTransport(theCascade);
  theFTFPanti->SetHighEnergyGenerator(theStringModel);

  //Ions
  //Binary Cascade + QMD + FTFP
  //theIonBC + QMD + theFTFPforIon
  G4ExcitationHandler* handler = new G4ExcitationHandler();
  G4PreCompoundModel* thePreCompound = new G4PreCompoundModel(handler);

  G4TheoFSGenerator* theFTFPforIonwtQMD;

  // Binary Cascade
  G4BinaryLightIonReaction* theIonBCwtQMD = new G4BinaryLightIonReaction(thePreCompound);
  theIonBCwtQMD->SetMinEnergy(0.0);
  theIonBCwtQMD->SetMaxEnergy(100*MeV);

  G4QMDReaction * theQMD= new G4QMDReaction;
  theQMD->SetMinEnergy(90*MeV);
  theQMD->SetMaxEnergy(10*GeV);

  // FTFP
  theFTFPforIonwtQMD = theFTFPwtBERT;
  theFTFPforIonwtQMD = new G4TheoFSGenerator("FTFP");

  theFTFPforIonwtQMD->SetMinEnergy( 9.*GeV );
  theFTFPforIonwtQMD->SetMaxEnergy( 100.*TeV );

/*
   theStringModel = new G4FTFModel;
   theStringDecay = new G4ExcitedStringDecay(theLund = new G4LundStringFragmentation);
   theStringModel->SetFragmentationModel(theStringDecay);

   theCascade = new G4GeneratorPrecompoundInterface;
   thePreEquilib = new G4PreCompoundModel(theHandler = new G4ExcitationHandler);
   theCascade->SetDeExcitation(thePreEquilib);

   theFTFPwtBERT->SetTransport(theCascade);
   theFTFPwtBERT->SetHighEnergyGenerator(theStringModel);

   theBertini = new G4CascadeInterface;
   theBertini->SetMinEnergy( 0.*GeV );
   theBertini->SetMaxEnergy( 6.*GeV );


    // this will be the model class for high energies
    G4TheoFSGenerator * theTheoModel = new G4TheoFSGenerator;


    // Pre equilibrium stage
    G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel;


    // a no-cascade generator-precompound interaface
    G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface;
            theCascade->SetDeExcitation(thePreEquilib);

    // here come the high energy parts
    // the string model; still not quite according to design - Explicite use of the forseen interfaces
    // will be tested and documented in this program by beta-02 at latest.
    G4VPartonStringModel * theStringModel;
    theStringModel = new G4QGSModel<G4QGSParticipants>;
    theTheoModel->SetTransport(theCascade);
    theTheoModel->SetHighEnergyGenerator(theStringModel);
    theTheoModel->SetMinEnergy(19*GeV);
    theTheoModel->SetMaxEnergy(100*TeV);

    G4BinaryCascade * theBC = new G4BinaryCascade;
    G4BinaryLightIonReaction * theIonBC= new G4BinaryLightIonReaction;
    theIonBC->SetMinEnergy(1*MeV);
    theIonBC->SetMaxEnergy(20*GeV);
*/


  G4TripathiCrossSection* TripathiCrossSection= new G4TripathiCrossSection;
  G4IonsShenCrossSection* aShen = new G4IonsShenCrossSection;
  // replace the default string fragmentation function (Lund) to QGSM
  //      G4VLongitudinalStringDecay * theFragmentation = new G4QGSMFragmentation;
  //      G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay(theFragmentation);
  //      theStringModel->SetFragmentationModel(theStringDecay);
  // done with the generator model (most of the above is also available as default)

  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess();
  theElasticProcess->RegisterMe(new G4HadronElastic());
  G4HadronElasticProcess* theElasticProcess1 = new G4HadronElasticProcess();

  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "pi+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4PionPlusInelasticProcess* theInelasticProcess =
                                new G4PionPlusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theBertini);
      theInelasticProcess->RegisterMe(theFTFPwtBERT);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "pi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4PionMinusInelasticProcess* theInelasticProcess =
                                new G4PionMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theBertini);
      theInelasticProcess->RegisterMe(theFTFPwtBERT);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonPlusInelasticProcess* theInelasticProcess =
                                  new G4KaonPlusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theBertini);
      theInelasticProcess->RegisterMe(theFTFPwtBERT);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon0S") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonZeroSInelasticProcess* theInelasticProcess =
                             new G4KaonZeroSInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theBertini);
      theInelasticProcess->RegisterMe(theFTFPwtBERT);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon0L") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonZeroLInelasticProcess* theInelasticProcess =
                             new G4KaonZeroLInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theBertini);
      theInelasticProcess->RegisterMe(theFTFPwtBERT);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonMinusInelasticProcess* theInelasticProcess =
                                 new G4KaonMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theBertini);
      theInelasticProcess->RegisterMe(theFTFPwtBERT);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "proton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4ProtonInelasticProcess* theInelasticProcess =
                                    new G4ProtonInelasticProcess("inelastic");
      //theInelasticProcess->RegisterMe(theBC);
      //theInelasticProcess->RegisterMe(theTheoModel);
      theInelasticProcess->RegisterMe(theBertini);
      theInelasticProcess->RegisterMe(theFTFPwtBERT);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_proton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiProtonInelasticProcess* theInelasticProcess =
                                new G4AntiProtonInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theFTFPanti);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "neutron") {
      // elastic scattering
      theElasticProcess1->RegisterMe(new G4HadronElastic());
      pmanager->AddDiscreteProcess(theElasticProcess1);
      // inelastic scattering
      G4NeutronInelasticProcess* theInelasticProcess =
                                    new G4NeutronInelasticProcess("inelastic");
      //theInelasticProcess->RegisterMe(theBC);
      //theInelasticProcess->RegisterMe(theTheoModel);
      theInelasticProcess->RegisterMe(theBertini);
      theInelasticProcess->RegisterMe(theFTFPwtBERT);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      // fission
      G4HadronFissionProcess* theFissionProcess =
                                    new G4HadronFissionProcess;
      G4LFission* theFissionModel = new G4LFission;
      theFissionProcess->RegisterMe(theFissionModel);
      pmanager->AddDiscreteProcess(theFissionProcess);
      // capture
      G4HadronCaptureProcess* theCaptureProcess =
                                    new G4HadronCaptureProcess;
      G4LCapture* theCaptureModel = new G4LCapture;
      theCaptureProcess->RegisterMe(theCaptureModel);
      pmanager->AddDiscreteProcess(theCaptureProcess);

    } else if (particleName == "anti_neutron") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiNeutronInelasticProcess* theInelasticProcess =
                               new G4AntiNeutronInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theFTFPanti);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "lambda") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LambdaInelasticProcess* theInelasticProcess =
                                    new G4LambdaInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theBertini);
      theInelasticProcess->RegisterMe(theFTFPwtBERT);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_lambda") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiLambdaInelasticProcess* theInelasticProcess =
                                new G4AntiLambdaInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theFTFPanti);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "sigma+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4SigmaPlusInelasticProcess* theInelasticProcess =
                                 new G4SigmaPlusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theBertini);
      theInelasticProcess->RegisterMe(theFTFPwtBERT);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "sigma-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4SigmaMinusInelasticProcess* theInelasticProcess =
                                 new G4SigmaMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theBertini);
      theInelasticProcess->RegisterMe(theFTFPwtBERT);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_sigma+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiSigmaPlusInelasticProcess* theInelasticProcess =
                             new G4AntiSigmaPlusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theFTFPanti);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_sigma-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiSigmaMinusInelasticProcess* theInelasticProcess =
                            new G4AntiSigmaMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theFTFPanti);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "xi0") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4XiZeroInelasticProcess* theInelasticProcess =
                            new G4XiZeroInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theBertini);
      theInelasticProcess->RegisterMe(theFTFPwtBERT);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "xi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4XiMinusInelasticProcess* theInelasticProcess =
                            new G4XiMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theBertini);
      theInelasticProcess->RegisterMe(theFTFPwtBERT);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "omega-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4OmegaMinusInelasticProcess* theInelasticProcess =
                            new G4OmegaMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theBertini);
      theInelasticProcess->RegisterMe(theFTFPwtBERT);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_xi0") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiXiZeroInelasticProcess* theInelasticProcess =
                            new G4AntiXiZeroInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theFTFPanti);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_xi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiXiMinusInelasticProcess* theInelasticProcess =
                            new G4AntiXiMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theFTFPanti);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_omega-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiOmegaMinusInelasticProcess* theInelasticProcess =
                            new G4AntiOmegaMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(theFTFPanti);
      pmanager->AddDiscreteProcess(theInelasticProcess);
 
    } else if (particleName == "deuteron") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4DeuteronInelasticProcess* theInelasticProcess =
                            new G4DeuteronInelasticProcess("inelastic");
      theInelasticProcess->AddDataSet(TripathiCrossSection);
      theInelasticProcess->AddDataSet(aShen);
      theInelasticProcess->RegisterMe(theIonBCwtQMD);
      theInelasticProcess->RegisterMe(theQMD);
      theInelasticProcess->RegisterMe(theFTFPforIonwtQMD);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "triton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4TritonInelasticProcess* theInelasticProcess =
                            new G4TritonInelasticProcess("inelastic");
      theInelasticProcess->AddDataSet(TripathiCrossSection);
      theInelasticProcess->AddDataSet(aShen);
      theInelasticProcess->RegisterMe(theIonBCwtQMD);
      theInelasticProcess->RegisterMe(theQMD);
      theInelasticProcess->RegisterMe(theFTFPforIonwtQMD);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "alpha") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AlphaInelasticProcess* theInelasticProcess =
                            new G4AlphaInelasticProcess("inelastic");
      theInelasticProcess->AddDataSet(TripathiCrossSection);
      theInelasticProcess->AddDataSet(aShen);
      theInelasticProcess->RegisterMe(theIonBCwtQMD);
      theInelasticProcess->RegisterMe(theQMD);
      theInelasticProcess->RegisterMe(theFTFPforIonwtQMD);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "GenericIon") {
      G4IonInelasticProcess* theInelasticProcess =
                            new G4IonInelasticProcess();
      theInelasticProcess->AddDataSet(TripathiCrossSection);
      theInelasticProcess->AddDataSet(aShen);
      //G4BinaryLightIonReaction * theGenIonBC= new G4BinaryLightIonReaction;
      //theGenIonBC->SetMinEnergy(0*MeV);
      //theGenIonBC->SetMaxEnergy(100*MeV);
      //G4QMDReaction * theQMD= new G4QMDReaction;
      //theQMD->SetMinEnergy(90*MeV);
      //theQMD->SetMaxEnergy(10*GeV);
      //theInelasticProcess->RegisterMe(theGenIonBC);
      //theInelasticProcess->RegisterMe(theQMD);
      theInelasticProcess->RegisterMe(theIonBCwtQMD);
      theInelasticProcess->RegisterMe(theQMD);
      theInelasticProcess->RegisterMe(theFTFPforIonwtQMD);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
  }
}

void Tst61PhysicsList::ConstructLeptHad()
{;}

#include "G4Decay.hh"
void Tst61PhysicsList::ConstructGeneral()
{
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) {
      pmanager ->AddProcess(theDecayProcess);
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}


void Tst61PhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "Tst61PhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue/mm << " (mm)" << G4endl;
  }
 //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
  SetCutsWithDefault();
}


