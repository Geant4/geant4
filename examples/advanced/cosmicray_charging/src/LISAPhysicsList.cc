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
// ********************************************************************
// *                                                                  *
// * cosmicray_charging advanced example for Geant4                   *
// * (adapted simulation of test-mass charging in the LISA mission)   *
// *                                                                  *
// * Henrique Araujo (h.araujo@imperial.ac.uk) & Peter Wass           *
// * Imperial College London                                          *
// *                                                                  *
// * LISAPhysicsList class                                            *
// *                                                                  *
// ********************************************************************
//
// HISTORY
// 22/02/2004: migrated from LISA-V04
//
// ********************************************************************

// Hadronics :
//           : PreCo for p, n at      0 < E < 70 MeV
//           : BiC   for p, n at 65 MeV < E < 6.1 GeV
//           : QGSP  for p, n at  6 GeV < E < 100 TeV
//           : BiC   for pi at        0 < E < 1.5 GeV
//           : LEP   for pi at  1.4 GeV < E < 6.1 GeV
//           : QGSP  for pi at    6 GeV < E < 100 TeV
//           : LEP for H2, H3, He4 at      0 < E < 100 MeV
//           : BiC for H2, H3, He4 at 80 MeV < E < 40 GeV
//           : BiC for He3, GenIon at      0 < E < 30 GeV



#include "LISAPhysicsList.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>


// Constructor /////////////////////////////////////////////////////////////
LISAPhysicsList::LISAPhysicsList() : G4VUserPhysicsList() {

  VerboseLevel = 1;
  SetVerboseLevel(VerboseLevel);
}


// Destructor //////////////////////////////////////////////////////////////
LISAPhysicsList::~LISAPhysicsList() 
{;}



////////////////////////////////////////////////////////////////////////////
// Construct Particles /////////////////////////////////////////////////////

#include "G4LeptonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

void LISAPhysicsList::ConstructParticle() {
  
  G4LeptonConstructor aC1;
  G4BaryonConstructor aC2;
  G4MesonConstructor aC3;
  G4BosonConstructor aC4;
  G4IonConstructor aC5;
  G4ShortLivedConstructor aC6;
  
  aC1.ConstructParticle();
  aC2.ConstructParticle();
  aC3.ConstructParticle();
  aC4.ConstructParticle();
  aC5.ConstructParticle();
  aC6.ConstructParticle();
  
}

/////////////////////////////////////////////////////////////////////////////
// Construct Processes //////////////////////////////////////////////////////

void LISAPhysicsList::ConstructProcess() {

  AddTransportation();

  ElectromagneticPhysics();

  HadronicPhysics();

  ElectroNuclearPhysics();

  GeneralPhysics();

}


/////////////////////////////////////////////////////////////////////////////
// Transportation ///////////////////////////////////////////////////////////

void LISAPhysicsList::AddTransportation() {

  G4VUserPhysicsList::AddTransportation();
  
}


/////////////////////////////////////////////////////////////////////////////
// Electromagnetic Processes ////////////////////////////////////////////////

// all charged particles
#include "G4MultipleScattering.hh"

// gamma
#include "G4LowEnergyRayleigh.hh" 
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyCompton.hh"  
#include "G4LowEnergyGammaConversion.hh" 

// e-
#include "G4LowEnergyIonisation.hh" 
#include "G4LowEnergyBremsstrahlung.hh" 

// e+
#include "G4eIonisation.hh" 
#include "G4eBremsstrahlung.hh" 
#include "G4eplusAnnihilation.hh"

// muons
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCaptureAtRest.hh"

// hadrons and ions
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"


void LISAPhysicsList::ElectromagneticPhysics() {


  G4cout << "Electromagnetic Physics" << G4endl;


   // processes

  G4LowEnergyPhotoElectric*  lowePhot = new G4LowEnergyPhotoElectric();
  G4LowEnergyIonisation*     loweIon  = new G4LowEnergyIonisation();
  G4LowEnergyBremsstrahlung* loweBrem = new G4LowEnergyBremsstrahlung();
  lowePhot->SetCutForLowEnSecPhotons(100*eV);
  loweIon ->SetCutForLowEnSecPhotons(100*eV);
  loweBrem->SetCutForLowEnSecPhotons(100*eV);

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName      = particle->GetParticleName();
    G4String particleType      = particle->GetParticleType();
    G4double particleCharge    = particle->GetPDGCharge();
    
    if (particleName == "gamma") {
      //gamma
      pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh());
      pmanager->AddDiscreteProcess(lowePhot);
      pmanager->AddDiscreteProcess(new G4LowEnergyCompton());
      pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion());

    } else if (particleName == "e-") {
      //electron
      G4MultipleScattering* aMultipleScattering = new G4MultipleScattering();
      // Modifying Facrange from default value (0.199) 
      // to improve backscattering fraction for electrons
      aMultipleScattering->SetFacrange(0.01);
      pmanager->AddProcess(aMultipleScattering,      -1, 1, 1);
      pmanager->AddProcess(loweIon,                  -1, 2, 2);
      pmanager->AddProcess(loweBrem,                 -1,-1, 3);

    } else if (particleName == "e+") {
      //positron
      G4MultipleScattering* aMultipleScattering = new G4MultipleScattering();
      pmanager->AddProcess(aMultipleScattering,      -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation(),      -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung(),  -1,-1, 3);
      pmanager->AddProcess(new G4eplusAnnihilation(), 0,-1, 4);      

    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon
      G4MultipleScattering* aMultipleScattering = new G4MultipleScattering();
      pmanager->AddProcess(aMultipleScattering,      -1, 1, 1);
      pmanager->AddProcess(new G4MuIonisation(),     -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung(), -1,-1, 3);
      pmanager->AddProcess(new G4MuPairProduction(), -1,-1, 4);
      if( particleName == "mu-" )
	pmanager->AddProcess(new G4MuonMinusCaptureAtRest(),0,-1,-1);
      
    } else if( particleName == "GenericIon" ) { 
      // ions
      pmanager->AddProcess(new G4MultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4ionIonisation,      -1, 2, 2);
      
    } else if (!particle->IsShortLived() &&
	       particleCharge != 0.0 && 
	       particleName   != "chargedgeantino") {
      // all other stable charged particles
      pmanager->AddProcess(new G4MultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,      -1, 2, 2);
    }
    
  }

}



///////////////////////////////////////////////////////////////////////////
// ElectroNuclear Physics /////////////////////////////////////////////////

// photonuclear and electronuclear reaction
#include "G4PhotoNuclearProcess.hh"
#include "G4ElectronNuclearProcess.hh"
#include "G4PositronNuclearProcess.hh"
#include "G4GammaNuclearReaction.hh"
#include "G4ElectroNuclearReaction.hh"

// CHIPS fragmentation model
#include "G4TheoFSGenerator.hh"
#include "G4StringChipsParticleLevelInterface.hh"
#include "G4QGSModel.hh"
#include "G4GammaParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

// muon photonuclear reaction
// #include "G4MuonNucleusProcess.hh"
// #include "G4MuNuclearInteraction.hh"


void LISAPhysicsList::ElectroNuclearPhysics() {

  G4cout << "ElectroNuclear Physics" << G4endl;

  // gamma
  G4PhotoNuclearProcess* thePhotoNuclearProcess = 
    new G4PhotoNuclearProcess("PhotoNuclear");
  // low energy
  G4GammaNuclearReaction* theGammaReaction = 
    new G4GammaNuclearReaction;
  theGammaReaction->SetMaxEnergy(3.*GeV);
  thePhotoNuclearProcess->RegisterMe(theGammaReaction);
  // high energy
  // CHIPS fragmentation for photonuclear reaction
  // Quark Gluon String model for high energy gammas
  G4TheoFSGenerator* theHEModel = new G4TheoFSGenerator;
  G4StringChipsParticleLevelInterface* theCascade = 
    new G4StringChipsParticleLevelInterface;
  G4QGSModel<G4GammaParticipants> theStringModel;
  theHEModel->SetTransport(theCascade);
  theHEModel->SetHighEnergyGenerator(&theStringModel);
  theHEModel->SetMinEnergy(3.*GeV);
  theHEModel->SetMaxEnergy(10*TeV);
  G4QGSMFragmentation theFragmentation;
  G4ExcitedStringDecay* theStringDecay = 
    new G4ExcitedStringDecay(&theFragmentation);
  theStringModel.SetFragmentationModel(theStringDecay);  
  thePhotoNuclearProcess->RegisterMe(theHEModel);
  //
  G4ProcessManager* pmanager = G4Gamma::Gamma()->GetProcessManager();
  pmanager->AddDiscreteProcess(thePhotoNuclearProcess);

  // e-
  pmanager = G4Electron::Electron()->GetProcessManager();
  G4ElectronNuclearProcess* theElectronNuclearProcess = 
    new G4ElectronNuclearProcess("ElectronNuclear");
  // see G4ElectroNuclearReaction.hh for defaults
  // 0 < E < 10 TeV
  G4ElectroNuclearReaction* theElectroReaction = 
    new G4ElectroNuclearReaction;
  theElectroReaction->SetMaxEnergy(10*TeV);
  theElectronNuclearProcess->RegisterMe(theElectroReaction);
  pmanager->AddDiscreteProcess(theElectronNuclearProcess);

  // e+
  pmanager = G4Positron::Positron()->GetProcessManager();
  G4PositronNuclearProcess* thePositronNuclearProcess = 
    new G4PositronNuclearProcess("PositronNuclear");
  // see G4ElectroNuclearReaction.hh for defaults
  // 0 < E < 10 TeV
  thePositronNuclearProcess->RegisterMe(theElectroReaction);
  pmanager->AddDiscreteProcess(thePositronNuclearProcess);


  // Muon spectra too soft for photonuclear reaction
  //    // mu-
  //    pmanager = G4MuonMinus::MuonMinus()->GetProcessManager();
  //    G4MuonNucleusProcess* theMuonNucleusProcess = 
  //      new G4MuonNucleusProcess("MuonNucleus");
  //    pmanager->AddDiscreteProcess(theMuonNucleusProcess);
  //    // mu+
  //    pmanager = G4MuonPlus::MuonPlus()->GetProcessManager();
  //    pmanager->AddDiscreteProcess(theMuonNucleusProcess);


}



/////////////////////////////////////////////////////////////////////////////
// Hadronic Physics /////////////////////////////////////////////////////////

#include "G4HadronElasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4HadronFissionProcess.hh"

// Inelastic Processes
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
// #include "G4IonInelasticProcess.hh"
#include "G4LambdaInelasticProcess.hh"
#include "G4AntiLambdaInelasticProcess.hh"
#include "G4OmegaMinusInelasticProcess.hh"
#include "G4AntiOmegaMinusInelasticProcess.hh"
#include "G4SigmaMinusInelasticProcess.hh"
#include "G4AntiSigmaMinusInelasticProcess.hh"
#include "G4SigmaPlusInelasticProcess.hh"
#include "G4AntiSigmaPlusInelasticProcess.hh"
#include "G4XiZeroInelasticProcess.hh"
#include "G4AntiXiZeroInelasticProcess.hh"
#include "G4XiMinusInelasticProcess.hh"
#include "G4AntiXiMinusInelasticProcess.hh"

// Low-energy Parameterised Models: 1 to 25 GeV
#include "G4LElastic.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"
// #include "G4LEProtonInelastic.hh"
#include "G4LEAntiProtonInelastic.hh"
// #include "G4LENeutronInelastic.hh"
#include "G4LEAntiNeutronInelastic.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"
#include "G4LELambdaInelastic.hh"
#include "G4LEAntiLambdaInelastic.hh"
#include "G4LEOmegaMinusInelastic.hh"
#include "G4LEAntiOmegaMinusInelastic.hh"
#include "G4LESigmaMinusInelastic.hh"
#include "G4LEAntiSigmaMinusInelastic.hh"
#include "G4LESigmaPlusInelastic.hh"
#include "G4LEAntiSigmaPlusInelastic.hh"
#include "G4LEXiZeroInelastic.hh"
#include "G4LEAntiXiZeroInelastic.hh"
#include "G4LEXiMinusInelastic.hh"
#include "G4LEAntiXiMinusInelastic.hh"
// neutrons
#include "G4LCapture.hh"
#include "G4LFission.hh"
// #include "G4ParaFissionModel.hh"

// High-energy Parameterised Models: 25 GeV to 10 TeV
//  #include "G4HEPionPlusInelastic.hh"
//  #include "G4HEPionMinusInelastic.hh"
//  #include "G4HEKaonPlusInelastic.hh"
//  #include "G4HEKaonZeroInelastic.hh"
//  #include "G4HEKaonZeroInelastic.hh"
//  #include "G4HEKaonMinusInelastic.hh"
//  #include "G4HEProtonInelastic.hh"
//  #include "G4HENeutronInelastic.hh"
#include "G4HEAntiProtonInelastic.hh"
#include "G4HEAntiNeutronInelastic.hh"
#include "G4HELambdaInelastic.hh"
#include "G4HEAntiLambdaInelastic.hh"
#include "G4HEOmegaMinusInelastic.hh"
#include "G4HEAntiOmegaMinusInelastic.hh"
#include "G4HESigmaMinusInelastic.hh"
#include "G4HEAntiSigmaMinusInelastic.hh"
#include "G4HESigmaPlusInelastic.hh"
#include "G4HEAntiSigmaPlusInelastic.hh"
#include "G4HEXiZeroInelastic.hh"
#include "G4HEAntiXiZeroInelastic.hh"
#include "G4HEXiMinusInelastic.hh"
#include "G4HEAntiXiMinusInelastic.hh"

// Neutron HP Models: Thermal to 19 MeV
// #include "G4NeutronHPElastic.hh"
// #include "G4NeutronHPElasticData.hh"
// #include "G4NeutronHPCapture.hh"
// #include "G4NeutronHPCaptureData.hh"
// #include "G4NeutronHPInelastic.hh"
// #include "G4NeutronHPInelasticData.hh"

// Stopping processes
#include "G4PiMinusAbsorptionAtRest.hh"
#include "G4KaonMinusAbsorptionAtRest.hh"
#include "G4AntiProtonAnnihilationAtRest.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"

// Generator models: HE
#include "G4TheoFSGenerator.hh"
#include "G4Evaporation.hh"
#include "G4CompetitiveFission.hh"
#include "G4FermiBreakUp.hh"
#include "G4StatMF.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"

// Kinetic Model
#include "G4BinaryCascade.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"


void LISAPhysicsList::HadronicPhysics() {


  G4cout << "Hadronic Physics" << G4endl;


  // **************************************************//
  // *** preparing inelastic reactions for hadrons *** // 
  // **************************************************//

  // high energy model for proton, neutron, pions and kaons
  G4TheoFSGenerator* theHEModel = new G4TheoFSGenerator;

  // all models for treatment of thermal nucleus 
  G4Evaporation* theEvaporation = new G4Evaporation;
  G4FermiBreakUp* theFermiBreakUp = new G4FermiBreakUp;
  G4StatMF* theMF = new G4StatMF;

  // evaporation logic
  G4ExcitationHandler* theHandler = new G4ExcitationHandler;
  theHandler->SetEvaporation(theEvaporation);
  theHandler->SetFermiModel(theFermiBreakUp);
  theHandler->SetMultiFragmentation(theMF);
  theHandler->SetMaxAandZForFermiBreakUp(12, 6);
  theHandler->SetMinEForMultiFrag(3.*MeV);

  // pre-equilibrium stage 
  G4PreCompoundModel* thePreEquilib = new G4PreCompoundModel(theHandler);
  thePreEquilib->SetMaxEnergy(70*MeV);

  // a no-cascade generator-precompound interaface
  G4GeneratorPrecompoundInterface* theCascade = 
    new G4GeneratorPrecompoundInterface;
  theCascade->SetDeExcitation(thePreEquilib);

  // QGSP model
  G4VPartonStringModel* theStringModel;
  theStringModel = new G4QGSModel<G4QGSParticipants>;
  theHEModel->SetTransport(theCascade);
  theHEModel->SetHighEnergyGenerator(theStringModel);
  theHEModel->SetMinEnergy(6*GeV);
  theHEModel->SetMaxEnergy(100*TeV);

  // Binary cascade for p, n
  G4BinaryCascade* theCasc = new G4BinaryCascade;
  theCasc->SetMinEnergy(65*MeV);
  theCasc->SetMaxEnergy(6.1*GeV);

  // fragmentation
  G4VLongitudinalStringDecay* theFragmentation = new G4QGSMFragmentation;
  G4ExcitedStringDecay* theStringDecay = 
    new G4ExcitedStringDecay(theFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);


  // pion inelastics
  // Binary cascade
  G4BinaryCascade* theCascForPi = new G4BinaryCascade;
  // theCascForPi->SetMinEnergy(65*MeV);
  theCascForPi->SetMinEnergy(0*MeV);
  theCascForPi->SetMaxEnergy(1.5*GeV);
  // LEP to fill the gap
  G4LEPionPlusInelastic* theLEPModelPiPlus = new G4LEPionPlusInelastic;
  theLEPModelPiPlus->SetMinEnergy(1.4*GeV);
  theLEPModelPiPlus->SetMaxEnergy(6.1*GeV);
  G4LEPionMinusInelastic* theLEPModelPiMinus = new G4LEPionMinusInelastic;
  theLEPModelPiMinus->SetMinEnergy(1.4*GeV);
  theLEPModelPiMinus->SetMaxEnergy(6.1*GeV);



  // *******************************************************//
  // *** preparing inelastic reactions for light nuclei *** // 
  // *******************************************************//

  // binary cascade for light nuclei
  // NOTE: Shen XS only up to 10 GeV/n;
  // pushing up to 40 GeV for alphas; will fail for d>20 GeV and t>30 GeV
  G4BinaryLightIonReaction * theIonBC= new G4BinaryLightIonReaction;
  theIonBC->SetMinEnergy(80*MeV);
  //  theIonBC->SetMaxEnergy(20*GeV);
  theIonBC->SetMaxEnergy(40*GeV);
  G4TripathiCrossSection* TripathiCrossSection= new G4TripathiCrossSection;
  G4IonsShenCrossSection* aShen = new G4IonsShenCrossSection;

  // deuteron
  G4LEDeuteronInelastic* theDIModel = new G4LEDeuteronInelastic;
  theDIModel->SetMaxEnergy(100*MeV);

  // triton
  G4LETritonInelastic* theTIModel = new G4LETritonInelastic;
  theTIModel->SetMaxEnergy(100*MeV);

  // alpha
  G4LEAlphaInelastic* theAIModel = new G4LEAlphaInelastic;
  theAIModel->SetMaxEnergy(100*MeV);

  // Generic Ion and He3
  // NOTE: Shen XS only up to 10 GeV/n;
  // pushing up to 30 GeV for He3;
  G4BinaryLightIonReaction* theGenIonBC= new G4BinaryLightIonReaction;
  theGenIonBC->SetMinEnergy(0*MeV);
  //  theGenIonBC->SetMaxEnergy(10*GeV);
  theGenIonBC->SetMaxEnergy(30*GeV);



  // *************************************//
  // *** preparing elastic scattering *** //
  // *************************************//

  // elastic process 
  // now producing nuclear recoils
  G4HadronElasticProcess* theElasticProcess = 
    new G4HadronElasticProcess("elastic");
  G4LElastic* theElasticModel = new G4LElastic;
  theElasticProcess->RegisterMe(theElasticModel);



  // *****************************************//
  // *** attaching processes to particles *** //
  // *****************************************//

  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "pi+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4PionPlusInelasticProcess* theInelasticProcess = 
	new G4PionPlusInelasticProcess("inelastic");
      // NOTE: PreCo crahes for Pi+
      // theInelasticProcess->RegisterMe(thePreEquilib);
      theInelasticProcess->RegisterMe(theCascForPi);
      theInelasticProcess->RegisterMe(theLEPModelPiPlus);
      theInelasticProcess->RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "pi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4PionMinusInelasticProcess* theInelasticProcess = 
	new G4PionMinusInelasticProcess("inelastic");
      // theInelasticProcess->RegisterMe(thePreEquilib);
      theInelasticProcess->RegisterMe(theCascForPi);
      theInelasticProcess->RegisterMe(theLEPModelPiMinus);
      theInelasticProcess->RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4PiMinusAbsorptionAtRest, ordDefault);

    } else if (particleName == "kaon+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonPlusInelasticProcess* theInelasticProcess = 
	new G4KaonPlusInelasticProcess("inelastic");
      G4LEKaonPlusInelastic* theLEInelasticModel = new G4LEKaonPlusInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      // G4HEKaonPlusInelastic* theHEInelasticModel = 
      // new G4HEKaonPlusInelastic;
      // theInelasticProcess->RegisterMe(theHEInelasticModel);
      theInelasticProcess->RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      
    } else if (particleName == "kaon0S") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonZeroSInelasticProcess* theInelasticProcess = 
	new G4KaonZeroSInelasticProcess("inelastic");
      G4LEKaonZeroSInelastic* theLEInelasticModel = new G4LEKaonZeroSInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      // G4HEKaonZeroInelastic* theHEInelasticModel = 
      // new G4HEKaonZeroInelastic;
      // theInelasticProcess->RegisterMe(theHEInelasticModel);
      theInelasticProcess->RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon0L") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonZeroLInelasticProcess* theInelasticProcess = 
	new G4KaonZeroLInelasticProcess("inelastic");
      G4LEKaonZeroLInelastic* theLEInelasticModel = new G4LEKaonZeroLInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      // G4HEKaonZeroInelastic* theHEInelasticModel = 
      // new G4HEKaonZeroInelastic;
      // theInelasticProcess->RegisterMe(theHEInelasticModel);
      theInelasticProcess->RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonMinusInelasticProcess* theInelasticProcess = 
	new G4KaonMinusInelasticProcess("inelastic");
      G4LEKaonMinusInelastic* theLEInelasticModel = new G4LEKaonMinusInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      // G4HEKaonMinusInelastic* theHEInelasticModel = 
      // new G4HEKaonMinusInelastic;
      // theInelasticProcess->RegisterMe(theHEInelasticModel);
      theInelasticProcess->RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4KaonMinusAbsorptionAtRest, ordDefault);

    } else if (particleName == "proton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4ProtonInelasticProcess* theInelasticProcess = 
	new G4ProtonInelasticProcess("inelastic");
      // G4LEProtonInelastic* theLEInelasticModel = new G4LEProtonInelastic;
      // theLEInelasticModel->SetMaxEnergy(25*GeV);
      // theInelasticProcess->RegisterMe(theLEInelasticModel);
      // G4HEProtonInelastic* theHEInelasticModel = new G4HEProtonInelastic;
      // theInelasticProcess->RegisterMe(theHEInelasticModel);
      theInelasticProcess->RegisterMe(thePreEquilib);
      theInelasticProcess->RegisterMe(theCasc);
      theInelasticProcess->RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_proton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiProtonInelasticProcess* theInelasticProcess = 
  	new G4AntiProtonInelasticProcess("inelastic");
      G4LEAntiProtonInelastic* theLEInelasticModel=new G4LEAntiProtonInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      G4HEAntiProtonInelastic* theHEInelasticModel = 
  	new G4HEAntiProtonInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4AntiProtonAnnihilationAtRest, ordDefault);

    } else if (particleName == "neutron") {
      // elastic scattering
      G4HadronElasticProcess* theNeutronElasticProcess = 
  	new G4HadronElasticProcess("elastic");
      G4LElastic* theElasticModel1 = new G4LElastic;
      //      theElasticModel1->SetMinEnergy(19*MeV);
      theNeutronElasticProcess->RegisterMe(theElasticModel1);
      // G4NeutronHPElastic* theElasticModel2 = new G4NeutronHPElastic;
      // theNeutronElasticProcess->RegisterMe(theElasticModel2);
      // G4CrossSectionDataStore* theElasticStore = 
      //   ((G4HadronElasticProcess*)theNeutronElasticProcess)
      //   ->GetCrossSectionDataStore();
      // G4NeutronHPElasticData* theElasticData = new G4NeutronHPElasticData;
      // theElasticStore->AddDataSet(theElasticData);
      pmanager->AddDiscreteProcess(theNeutronElasticProcess);

      // inelastic scattering
      G4NeutronInelasticProcess* theInelasticProcess =
  	new G4NeutronInelasticProcess("inelastic");
      //       // inelastic model 1: HP
      //       G4NeutronHPInelastic* theInelasticModel1 =
      //    	new G4NeutronHPInelastic;
      //       theInelasticProcess->RegisterMe(theInelasticModel1);
      //       G4CrossSectionDataStore* theInelasticStore = 
      //    	((G4HadronInelasticProcess*)theInelasticProcess)
      //    	->GetCrossSectionDataStore();
      //       G4NeutronHPInelasticData* theInelasticData1 = 
      //    	new G4NeutronHPInelasticData;
      //       theInelasticStore->AddDataSet(theInelasticData1);
      // inelastic model 2: HK
      // G4LENeutronInelastic* theInelasticModel2 = new G4LENeutronInelastic;
      // theInelasticModel2->SetMinEnergy(19*MeV);
      // theInelasticModel2->SetMaxEnergy(25*GeV);
      // theInelasticProcess->RegisterMe(theInelasticModel2);
      // using pre-equilibrium instead of neutronHP
      theInelasticProcess->RegisterMe(thePreEquilib);
      theInelasticProcess->RegisterMe(theCasc);
      theInelasticProcess->RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

      // capture
      G4HadronCaptureProcess* theCaptureProcess =
    	new G4HadronCaptureProcess("capture");
      G4LCapture* theCaptureModel = new G4LCapture;
      // theCaptureModel->SetMinEnergy(19*MeV);
      theCaptureProcess->RegisterMe(theCaptureModel);
      // G4NeutronHPCapture* theLENeutronCaptureModel = new G4NeutronHPCapture;
      // theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
      // G4CrossSectionDataStore* theCaptureStore = 
      // ((G4HadronCaptureProcess*)theCaptureProcess)->
      // GetCrossSectionDataStore();
      // G4NeutronHPCaptureData* theCaptureData = new G4NeutronHPCaptureData;
      // theCaptureStore->AddDataSet(theCaptureData);
      pmanager->AddDiscreteProcess(theCaptureProcess);

      // fission
      G4HadronFissionProcess* theFissionProcess =
	new G4HadronFissionProcess;
      G4LFission* theFissionModel = new G4LFission;
      // G4ParaFissionModel* theFissionModel = new G4ParaFissionModel;
      theFissionProcess->RegisterMe(theFissionModel);
      pmanager->AddDiscreteProcess(theFissionProcess);

    } else if (particleName == "anti_neutron") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiNeutronInelasticProcess* theInelasticProcess = 
  	new G4AntiNeutronInelasticProcess("inelastic");
      G4LEAntiNeutronInelastic* theLEInelasticModel = 
  	new G4LEAntiNeutronInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      G4HEAntiNeutronInelastic* theHEInelasticModel = 
  	new G4HEAntiNeutronInelastic;
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4AntiNeutronAnnihilationAtRest,ordDefault);

    } else if (particleName == "deuteron") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4DeuteronInelasticProcess* theIPdeuteron = 
	new G4DeuteronInelasticProcess;
      theIPdeuteron->AddDataSet(TripathiCrossSection);
      theIPdeuteron->AddDataSet(aShen);
      theIPdeuteron->RegisterMe(theDIModel);
      theIPdeuteron->RegisterMe(theIonBC);
      pmanager->AddDiscreteProcess(theIPdeuteron);

    } else if (particleName == "triton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4TritonInelasticProcess* theIPtriton = 
	new G4TritonInelasticProcess;
      theIPtriton->AddDataSet(TripathiCrossSection);
      theIPtriton->AddDataSet(aShen);
      theIPtriton->RegisterMe(theTIModel);
      theIPtriton->RegisterMe(theIonBC);
      pmanager->AddDiscreteProcess(theIPtriton);

    } else if (particleName == "alpha") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AlphaInelasticProcess* theIPalpha = 
	new G4AlphaInelasticProcess;
      theIPalpha->AddDataSet(TripathiCrossSection);
      theIPalpha->AddDataSet(aShen);
      theIPalpha->RegisterMe(theAIModel);
      theIPalpha->RegisterMe(theIonBC);
      pmanager->AddDiscreteProcess(theIPalpha);

    } else if (particleName == "He3") {
      // NOTE elastic scattering does not stick to He3!
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theIPHe3 =
	new G4HadronInelasticProcess("He3Inelastic", G4He3::He3());
      theIPHe3->AddDataSet(TripathiCrossSection);
      theIPHe3->AddDataSet(aShen);
      theIPHe3->RegisterMe(theGenIonBC);
      pmanager->AddDiscreteProcess(theIPHe3);

    } else if (particleName == "GenericIon") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theIPGenericIon =
	new G4HadronInelasticProcess("IonInel",G4GenericIon::GenericIon());
      theIPGenericIon->AddDataSet(TripathiCrossSection);
      theIPGenericIon->AddDataSet(aShen);
      theIPGenericIon->RegisterMe(theGenIonBC);
      pmanager->AddDiscreteProcess(theIPGenericIon);

    } else if (particleName == "lambda") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LambdaInelasticProcess* theInelasticProcess = 
	new G4LambdaInelasticProcess("inelastic");
      G4LELambdaInelastic* theLEInelasticModel = new G4LELambdaInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      G4HELambdaInelastic* theHEInelasticModel = 
	new G4HELambdaInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_lambda") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiLambdaInelasticProcess* theInelasticProcess = 
	new G4AntiLambdaInelasticProcess("inelastic");
      G4LEAntiLambdaInelastic* theLEInelasticModel=new G4LEAntiLambdaInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      G4HEAntiLambdaInelastic* theHEInelasticModel = 
	new G4HEAntiLambdaInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "omega-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4OmegaMinusInelasticProcess* theInelasticProcess = 
	new G4OmegaMinusInelasticProcess("inelastic");
      G4LEOmegaMinusInelastic* theLEInelasticModel = 
	new G4LEOmegaMinusInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      G4HEOmegaMinusInelastic* theHEInelasticModel = 
	new G4HEOmegaMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_omega-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiOmegaMinusInelasticProcess* theInelasticProcess = 
	new G4AntiOmegaMinusInelasticProcess("inelastic");
      G4LEAntiOmegaMinusInelastic* theLEInelasticModel = 
	new G4LEAntiOmegaMinusInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      G4HEAntiOmegaMinusInelastic* theHEInelasticModel = 
	new G4HEAntiOmegaMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "sigma-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4SigmaMinusInelasticProcess* theInelasticProcess = 
	new G4SigmaMinusInelasticProcess("inelastic");
      G4LESigmaMinusInelastic* theLEInelasticModel = 
	new G4LESigmaMinusInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      G4HESigmaMinusInelastic* theHEInelasticModel = 
	new G4HESigmaMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_sigma-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiSigmaMinusInelasticProcess* theInelasticProcess = 
	new G4AntiSigmaMinusInelasticProcess("inelastic");
      G4LEAntiSigmaMinusInelastic* theLEInelasticModel = 
	new G4LEAntiSigmaMinusInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      G4HEAntiSigmaMinusInelastic* theHEInelasticModel = 
	new G4HEAntiSigmaMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "sigma+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4SigmaPlusInelasticProcess* theInelasticProcess = 
	new G4SigmaPlusInelasticProcess("inelastic");
      G4LESigmaPlusInelastic* theLEInelasticModel = 
	new G4LESigmaPlusInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      G4HESigmaPlusInelastic* theHEInelasticModel = 
	new G4HESigmaPlusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_sigma+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiSigmaPlusInelasticProcess* theInelasticProcess = 
	new G4AntiSigmaPlusInelasticProcess("inelastic");
      G4LEAntiSigmaPlusInelastic* theLEInelasticModel = 
	new G4LEAntiSigmaPlusInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      G4HEAntiSigmaPlusInelastic* theHEInelasticModel = 
	new G4HEAntiSigmaPlusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "xi0") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4XiZeroInelasticProcess* theInelasticProcess = 
	new G4XiZeroInelasticProcess("inelastic");
      G4LEXiZeroInelastic* theLEInelasticModel = 
	new G4LEXiZeroInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      G4HEXiZeroInelastic* theHEInelasticModel = 
	new G4HEXiZeroInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_xi0") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiXiZeroInelasticProcess* theInelasticProcess = 
	new G4AntiXiZeroInelasticProcess("inelastic");
      G4LEAntiXiZeroInelastic* theLEInelasticModel = 
	new G4LEAntiXiZeroInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      G4HEAntiXiZeroInelastic* theHEInelasticModel = 
	new G4HEAntiXiZeroInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "xi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4XiMinusInelasticProcess* theInelasticProcess = 
	new G4XiMinusInelasticProcess("inelastic");
      G4LEXiMinusInelastic* theLEInelasticModel = 
	new G4LEXiMinusInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      G4HEXiMinusInelastic* theHEInelasticModel = 
	new G4HEXiMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_xi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiXiMinusInelasticProcess* theInelasticProcess = 
	new G4AntiXiMinusInelasticProcess("inelastic");
      G4LEAntiXiMinusInelastic* theLEInelasticModel = 
	new G4LEAntiXiMinusInelastic;
      theLEInelasticModel->SetMaxEnergy(25*GeV);
      theInelasticProcess->RegisterMe(theLEInelasticModel);
      G4HEAntiXiMinusInelastic* theHEInelasticModel = 
	new G4HEAntiXiMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
  }
}


// Decays ///////////////////////////////////////////////////////////////////
#include "G4Decay.hh"
// #include "G4RadioactiveDecay.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"

void LISAPhysicsList::GeneralPhysics() {

  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay("decay");
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



// Cuts /////////////////////////////////////////////////////////////////////
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

void LISAPhysicsList::SetCuts() {

  // low energy limit
  G4double lowlimit=250*eV;
  G4ProductionCutsTable::GetProductionCutsTable()
    ->SetEnergyRange(lowlimit, 100.*GeV);

  // default cuts for world volume
  G4double cutValue = 2.0*mm;
  SetCutValue(cutValue,"gamma");
  SetCutValue(cutValue,"e-");
  SetCutValue(cutValue,"e+");

  // cuts per region: inertial sensor (250 eV for gamma, e-, e+)
  G4Region* region = G4RegionStore::GetInstance()->GetRegion("sensor");
  G4ProductionCuts* cuts = new G4ProductionCuts();
  cuts->SetProductionCut(0.050*micrometer);
  region->SetProductionCuts(cuts);
  
  if (verboseLevel>0) DumpCutValuesTable();

}

