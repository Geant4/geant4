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
// $Id: Tst38PhysicsList.cc,v 1.1 2006-05-02 10:27:31 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "Tst38PhysicsList.hh"

//Tst38PhysicsList::Tst38PhysicsList():  G4VModularPhysicsList()
Tst38PhysicsList::Tst38PhysicsList():  G4VUserPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(1);

  // General Physics
  //RegisterPhysics( new Tst38GeneralPhysics("general") );

  // EM Physics
  //RegisterPhysics( new Tst38EMPhysics("standard EM"));

  // Muon Physics
  //RegisterPhysics(  new Tst38MuonPhysics("muon"));

  // Hadron Physics
  //RegisterPhysics(  new Tst38HadronPhysics("hadron"));

  // Ion Physics
  //RegisterPhysics( new Tst38IonPhysics("ion"));

}

Tst38PhysicsList::~Tst38PhysicsList() {}

void Tst38PhysicsList::SetCuts() {SetCutsWithDefault();}


void Tst38PhysicsList::ConstructParticle()
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

void Tst38PhysicsList::ConstructAllBosons()
{
  // Construct all bosons
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst38PhysicsList::ConstructAllLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst38PhysicsList::ConstructAllMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst38PhysicsList::ConstructAllBaryons()
{
  //  Construct all barions
  G4BaryonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst38PhysicsList::ConstructAllIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void Tst38PhysicsList::ConstructAllShortLiveds()
{
  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void Tst38PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructLeptHad();
  ConstructHad();
  ConstructGeneral();
}

void Tst38PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") 
    {
      // Construct processes for gamma
      // ........The OLD (GHAR) usage of the CHIPS Photo-nuclear reactions .............
      //G4PhotoNuclearProcess* thePhotoNuclearProcess = new G4PhotoNuclearProcess;
      //G4QGSMFragmentation* theFragmentation = new G4QGSMFragmentation;
						//G4QGSModel<G4GammaParticipants>*theStringModel=new G4QGSModel<G4GammaParticipants>;
      //G4GammaNuclearReaction* theGammaReaction = new G4GammaNuclearReaction;
      //G4TheoFSGenerator* theModel = new G4TheoFSGenerator;
      //G4StringChipsParticleLevelInterface* theCascade =
      //                                           new G4StringChipsParticleLevelInterface;
						//theModel->SetTransport(theCascade);
      //theModel->SetHighEnergyGenerator(theStringModel);
      //G4ExcitedStringDecay* theStringDecay = new G4ExcitedStringDecay(theFragmentation);
      //theStringModel->SetFragmentationModel(theStringDecay);
      //theGammaReaction->SetMaxEnergy(3.5*GeV);
      //thePhotoNuclearProcess->RegisterMe(theGammaReaction);
      //theModel->SetMinEnergy(3.*GeV);
      //theModel->SetMaxEnergy(100*TeV);
      //thePhotoNuclearProcess->RegisterMe(theModel);
      //pmanager->AddDiscreteProcess(thePhotoNuclearProcess);
      // ........The NEW (native) usage of the CHIPS Photo-nuclear reactions .............
      pmanager->AddDiscreteProcess(new G4QCollision);

      // Pure Electromagnetic Processes
      pmanager->AddDiscreteProcess(new G4GammaConversion());     // Pair Production
      pmanager->AddDiscreteProcess(new G4ComptonScattering());   // Compton Effect     
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect()); // Photo Effect

    }
    else if (particleName == "e-")
    {
      // Construct processes for electron
      // ........The OLD (GHAR) usage of the CHIPS Photo-nuclear reactions .............
      //G4ElectronNuclearProcess* theElectronNuclearProcess = new G4ElectronNuclearProcess;
      //G4ElectroNuclearReaction* theElectronReaction = new G4ElectroNuclearReaction;
      //theElectronNuclearProcess->RegisterMe(theElectronReaction);
      //theElectronNuclearProcess->BiasCrossSectionByFactor(1000);
      //pmanager->AddDiscreteProcess(theElectronNuclearProcess);
      // ........The NEW (native) usage of the CHIPS electro-nuclear reactions ............
      pmanager->AddDiscreteProcess(new G4QCollision);

      // Pure Electromagnetic Processes
      pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);// ElectronMultipleScattering
      pmanager->AddProcess(new G4eIonisation(),-1,2,2);       // Electron Ionisation							
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);  // Electron BremsStrahlung   
  
    }
    else if (particleName == "e+")
    {
      // Construct processes for positron
      // ........The OLD (GHAR) usage of the CHIPS Photo-nuclear reactions .............
      //G4PositronNuclearProcess* thePositronNuclearProcess = new G4PositronNuclearProcess;
      //G4ElectroNuclearReaction* thePositronReaction = new G4ElectroNuclearReaction;
      //thePositronNuclearProcess->RegisterMe(thePositronReaction);
      //thePositronNuclearProcess->BiasCrossSectionByFactor(1000);
      //pmanager->AddDiscreteProcess(thePositronNuclearProcess);
      // ........The NEW (native) usage of the CHIPS positron-nuclear reactions ............
      pmanager->AddDiscreteProcess(new G4QCollision);

      // Pure Electromagnetic Processes
      pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);// PositronMultipleScattering
      pmanager->AddProcess(new G4eIonisation(),-1,2,2);							// Positron Ionisation
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);  // Positron BremsStrahlung
      G4eplusAnnihilation* theAnnihilation = new G4eplusAnnihilation;
      pmanager->AddDiscreteProcess(theAnnihilation);     // Positron Annihilation on Flight
      pmanager->AddRestProcess(theAnnihilation);         // Positron Annihilation at Rest
  
    }
    else if( particleName == "mu+" || particleName == "mu-"    )
    {
      // Construct processes for muon+/-
      pmanager->AddProcess(new G4MultipleScattering(),-1,1,1); // Mu Multiple Scattering
      pmanager->AddProcess(new G4MuIonisation(),-1,2,2);       // Mu Ionization
      pmanager->AddProcess(new G4MuBremsstrahlung(),-1,-1,3);  // Mu Bremsstrahlung
      pmanager->AddProcess(new G4MuPairProduction(),-1,-1,4);  // Mu Pair Production      
      // ........The NEW (native) usage of the CHIPS Muo-nuclear reactions .............
      pmanager->AddDiscreteProcess(new G4QCollision);
      //@@ the G4MuonMinusCaptureAtRest process is switched off
     
    }
    else if( particleName == "tau+" || particleName == "tau-"    )
    {
      // Construct processes for muon+/-
      pmanager->AddProcess(new G4MultipleScattering(),-1,1,1); // Tau Multiple Scattering

      // ........The NEW (native) usage of the CHIPS Muo-nuclear reactions .............
      pmanager->AddDiscreteProcess(new G4QCollision);
      //@@ the G4MuonMinusCaptureAtRest process for tau is switched off

    }
    else if( particleName == "GenericIon" )
    {
      // Construct processes for ions
      pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);
      pmanager->AddProcess(new G4hIonisation(),-1,2,2); 
    }
    else
    { 
      if (particle->GetPDGCharge() && (particle->GetParticleName() != "chargedgeantino") &&
                                                                !particle->IsShortLived() )
      {  
        // short lived particles except geantino
        pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);
        pmanager->AddProcess(new G4hIonisation(),-1,2,2);       
      }
    }
  }
}
//
// Makes discrete physics processes for the hadrons, 
//

void Tst38PhysicsList::ConstructHad()
{
  // this will be the model class for high energies
  G4TheoFSGenerator * theTheoModel = new G4TheoFSGenerator;
         
  // a cascade interface interface to chiral invriant phase space decay
  // at particle level.
  G4StringChipsParticleLevelInterface* theCascade =
                                                 new G4StringChipsParticleLevelInterface();
	
  // here come the high energy parts
  // the string model; still not quite according to design - Explicite use of the forseen
  // interfaces  will be tested and documented in this program by beta-02 at latest.
  G4VPartonStringModel* theStringModel;
  theStringModel = new G4QGSModel<G4QGSParticipants>;
  theTheoModel->SetTransport(theCascade);
  theTheoModel->SetHighEnergyGenerator(theStringModel);
  theTheoModel->SetMinEnergy(19*GeV);
  theTheoModel->SetMaxEnergy(100*TeV);

  G4VLongitudinalStringDecay * theFragmentation = new G4LundStringFragmentation;
  G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay(theFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);

  // done with the generator model (most of the above is also available as default)
  G4HadronElasticProcess* theElasticProcess =  new G4HadronElasticProcess; //@@
  G4LElastic* theElasticModel = new G4LElastic;                            //@@
  theElasticProcess->RegisterMe(theElasticModel);
  G4HadronElasticProcess* theElasticProcess1 =  new G4HadronElasticProcess;
  theParticleIterator->reset();
  while ((*theParticleIterator)())
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "pi+")
    {
       pmanager->AddDiscreteProcess(theElasticProcess);
       G4PionPlusInelasticProcess* theInelasticProcess = 
                                new G4PionPlusInelasticProcess("inelastic");
       G4LEPionPlusInelastic* theInelasticModel =  new G4LEPionPlusInelastic;
       theInelasticProcess->RegisterMe(theInelasticModel);
       theInelasticProcess->RegisterMe(theTheoModel);
       pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "pi-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4PionMinusInelasticProcess* theInelasticProcess = 
                                new G4PionMinusInelasticProcess("inelastic");
      G4LEPionMinusInelastic* theInelasticModel =  new G4LEPionMinusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "kaon+")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonPlusInelasticProcess* theInelasticProcess = 
                                  new G4KaonPlusInelasticProcess("inelastic");
      G4LEKaonPlusInelastic* theInelasticModel = new G4LEKaonPlusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      new G4HEKaonPlusInelastic; // @@
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "kaon0S")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonZeroSInelasticProcess* theInelasticProcess = 
                             new G4KaonZeroSInelasticProcess("inelastic");
      G4LEKaonZeroSInelastic* theInelasticModel =  new G4LEKaonZeroSInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "kaon0L")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonZeroLInelasticProcess* theInelasticProcess = 
                             new G4KaonZeroLInelasticProcess("inelastic");
      G4LEKaonZeroLInelastic* theInelasticModel =  new G4LEKaonZeroLInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      new G4HEKaonZeroInelastic;
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "kaon-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonMinusInelasticProcess* theInelasticProcess = 
                                 new G4KaonMinusInelasticProcess("inelastic");
      G4LEKaonMinusInelastic* theInelasticModel =  new G4LEKaonMinusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "proton")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4ProtonInelasticProcess* theInelasticProcess = 
                                    new G4ProtonInelasticProcess("inelastic");
      G4LEProtonInelastic* theInelasticModel = new G4LEProtonInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_proton")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiProtonInelasticProcess* theInelasticProcess = 
                                new G4AntiProtonInelasticProcess("inelastic");
      G4LEAntiProtonInelastic* theInelasticModel =  new G4LEAntiProtonInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEAntiProtonInelastic* theHEInelasticModel =  new G4HEAntiProtonInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "neutron")
    {     
      // elastic scattering
      G4LElastic* theElasticModel1 = new G4LElastic;
      theElasticProcess1->RegisterMe(theElasticModel1);
      pmanager->AddDiscreteProcess(theElasticProcess1);
      // inelastic scattering
      G4NeutronInelasticProcess* theInelasticProcess = 
                                    new G4NeutronInelasticProcess("inelastic");
      G4LENeutronInelastic* theInelasticModel = new G4LENeutronInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      // fission
      G4HadronFissionProcess* theFissionProcess = new G4HadronFissionProcess;
      G4LFission* theFissionModel = new G4LFission;
      theFissionProcess->RegisterMe(theFissionModel);
      pmanager->AddDiscreteProcess(theFissionProcess);
      // capture
      G4HadronCaptureProcess* theCaptureProcess = new G4HadronCaptureProcess;
      G4LCapture* theCaptureModel = new G4LCapture;
      theCaptureProcess->RegisterMe(theCaptureModel);
      pmanager->AddDiscreteProcess(theCaptureProcess);
    }  
    else if (particleName == "anti_neutron")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiNeutronInelasticProcess* theInelasticProcess = 
                               new G4AntiNeutronInelasticProcess("inelastic");
      G4LEAntiNeutronInelastic* theInelasticModel =  new G4LEAntiNeutronInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEAntiNeutronInelastic* theHEInelasticModel =  new G4HEAntiNeutronInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "lambda")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LambdaInelasticProcess* theInelasticProcess = 
                                    new G4LambdaInelasticProcess("inelastic");
      G4LELambdaInelastic* theInelasticModel = new G4LELambdaInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HELambdaInelastic* theHEInelasticModel = new G4HELambdaInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_lambda")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiLambdaInelasticProcess* theInelasticProcess = 
                                new G4AntiLambdaInelasticProcess("inelastic");
      G4LEAntiLambdaInelastic* theInelasticModel =  new G4LEAntiLambdaInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEAntiLambdaInelastic* theHEInelasticModel =  new G4HEAntiLambdaInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "sigma+")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4SigmaPlusInelasticProcess* theInelasticProcess = 
                                 new G4SigmaPlusInelasticProcess("inelastic");
      G4LESigmaPlusInelastic* theInelasticModel =  new G4LESigmaPlusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HESigmaPlusInelastic* theHEInelasticModel =  new G4HESigmaPlusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "sigma-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4SigmaMinusInelasticProcess* theInelasticProcess = 
                                 new G4SigmaMinusInelasticProcess("inelastic");
      G4LESigmaMinusInelastic* theInelasticModel =  new G4LESigmaMinusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HESigmaMinusInelastic* theHEInelasticModel =  new G4HESigmaMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_sigma+")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiSigmaPlusInelasticProcess* theInelasticProcess = 
                             new G4AntiSigmaPlusInelasticProcess("inelastic");
      G4LEAntiSigmaPlusInelastic* theInelasticModel =  new G4LEAntiSigmaPlusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEAntiSigmaPlusInelastic* theHEInelasticModel =  new G4HEAntiSigmaPlusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_sigma-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiSigmaMinusInelasticProcess* theInelasticProcess = 
                            new G4AntiSigmaMinusInelasticProcess("inelastic");
      G4LEAntiSigmaMinusInelastic* theInelasticModel =  new G4LEAntiSigmaMinusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEAntiSigmaMinusInelastic* theHEInelasticModel =  new G4HEAntiSigmaMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "xi0")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4XiZeroInelasticProcess* theInelasticProcess = 
                            new G4XiZeroInelasticProcess("inelastic");
      G4LEXiZeroInelastic* theInelasticModel =  new G4LEXiZeroInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEXiZeroInelastic* theHEInelasticModel =  new G4HEXiZeroInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "xi-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4XiMinusInelasticProcess* theInelasticProcess = 
                            new G4XiMinusInelasticProcess("inelastic");
      G4LEXiMinusInelastic* theInelasticModel =  new G4LEXiMinusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEXiMinusInelastic* theHEInelasticModel =  new G4HEXiMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_xi0")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiXiZeroInelasticProcess* theInelasticProcess = 
                            new G4AntiXiZeroInelasticProcess("inelastic");
      G4LEAntiXiZeroInelastic* theInelasticModel =  new G4LEAntiXiZeroInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEAntiXiZeroInelastic* theHEInelasticModel =  new G4HEAntiXiZeroInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_xi-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiXiMinusInelasticProcess* theInelasticProcess = 
                            new G4AntiXiMinusInelasticProcess("inelastic");
      G4LEAntiXiMinusInelastic* theInelasticModel =  new G4LEAntiXiMinusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEAntiXiMinusInelastic* theHEInelasticModel =  new G4HEAntiXiMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "deuteron")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4DeuteronInelasticProcess* theInelasticProcess = 
                            new G4DeuteronInelasticProcess("inelastic");
      G4LEDeuteronInelastic* theInelasticModel =  new G4LEDeuteronInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "triton")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4TritonInelasticProcess* theInelasticProcess = 
                            new G4TritonInelasticProcess("inelastic");
      G4LETritonInelastic* theInelasticModel =  new G4LETritonInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "alpha")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AlphaInelasticProcess* theInelasticProcess = 
                            new G4AlphaInelasticProcess("inelastic");
      G4LEAlphaInelastic* theInelasticModel =  new G4LEAlphaInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "omega-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4OmegaMinusInelasticProcess* theInelasticProcess = 
                            new G4OmegaMinusInelasticProcess("inelastic");
      G4LEOmegaMinusInelastic* theInelasticModel =  new G4LEOmegaMinusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEOmegaMinusInelastic* theHEInelasticModel =  new G4HEOmegaMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_omega-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiOmegaMinusInelasticProcess* theInelasticProcess = 
                            new G4AntiOmegaMinusInelasticProcess("inelastic");
      G4LEAntiOmegaMinusInelastic* theInelasticModel =  new G4LEAntiOmegaMinusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEAntiOmegaMinusInelastic* theHEInelasticModel =  new G4HEAntiOmegaMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
  }
}

void Tst38PhysicsList::ConstructLeptHad() {;}                      //@@ Absolete

void Tst38PhysicsList::ConstructGeneral()
{
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle))
    { 
      pmanager ->AddProcess(theDecayProcess);
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}
