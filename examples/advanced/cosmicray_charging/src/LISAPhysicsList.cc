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
// 09/08/2004: Removed call by pointer of hadronics classes
// 09/08/2004: Added MuNuclear interaction
// 08/12/2005: changed particle construction
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

  G4LeptonConstructor lepton;
  lepton.ConstructParticle();
  
  G4BosonConstructor boson;
  boson.ConstructParticle();
  
  G4MesonConstructor meson;
  meson.ConstructParticle();
  
  G4BaryonConstructor baryon;
  baryon.ConstructParticle();
  
  G4ShortLivedConstructor shortLived;
  shortLived.ConstructParticle();
  
  G4IonConstructor ion;
  ion.ConstructParticle();

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
void LISAPhysicsList::ElectroNuclearPhysics() {

  G4cout << "ElectroNuclear Physics" << G4endl;

  // gamma
  G4ProcessManager* pmanager = G4Gamma::Gamma()->GetProcessManager();
  // low energy
  theGammaReaction = new G4GammaNuclearReaction;
  theGammaReaction->SetMaxEnergy(3.5*GeV);
  thePhotoNuclearProcess.RegisterMe(theGammaReaction);
  // high energy
  theHEModel_PN = new G4TheoFSGenerator;
  theCascade_PN = new G4StringChipsParticleLevelInterface;
  theHEModel_PN->SetTransport(theCascade_PN);
  theHEModel_PN->SetHighEnergyGenerator(theStringModel_PN);
  theStringDecay_PN = new G4ExcitedStringDecay(&theFragmentation_PN);

  theStringModel_PN = new G4QGSModel<G4GammaParticipants>;
  theStringModel_PN -> SetFragmentationModel(theStringDecay_PN);
  theHEModel_PN->SetMinEnergy(3.*GeV);
  theHEModel_PN->SetMaxEnergy(100*TeV);
  thePhotoNuclearProcess.RegisterMe(theHEModel_PN);
  pmanager->AddDiscreteProcess(&thePhotoNuclearProcess);
  
  // e-
  pmanager = G4Electron::Electron()->GetProcessManager();
  // see G4ElectroNuclearReaction.hh for defaults
  // 0 < E < 10 TeV
  theElectroReaction = new G4ElectroNuclearReaction;
  theElectroReaction->SetMaxEnergy(10*TeV);
  theElectronNuclearProcess.RegisterMe(theElectroReaction);
  pmanager->AddDiscreteProcess(&theElectronNuclearProcess);

  // e+
  pmanager = G4Positron::Positron()->GetProcessManager();
  // see G4ElectroNuclearReaction.hh for defaults
  // 0 < E < 10 TeV
  thePositronNuclearProcess.RegisterMe(theElectroReaction);
  pmanager->AddDiscreteProcess(&thePositronNuclearProcess);

  // Mu-nuclear reaction
  // mu-
  pmanager = G4MuonMinus::MuonMinus()->GetProcessManager();
  pmanager->AddDiscreteProcess(&theMuMinusNuclearInteraction);
  // mu+
  pmanager = G4MuonPlus::MuonPlus()->GetProcessManager();
  pmanager->AddDiscreteProcess(&theMuPlusNuclearInteraction);


}



/////////////////////////////////////////////////////////////////////////////
// Hadronic Physics /////////////////////////////////////////////////////////
void LISAPhysicsList::HadronicPhysics() {


  G4cout << "Hadronic Physics" << G4endl;


  // **************************************************//
  // *** preparing inelastic reactions for hadrons *** // 
  // **************************************************//
  //
  // high energy model for proton, neutron, pions and kaons
  theHEModel = new G4TheoFSGenerator;
  // all models for treatment of thermal nucleus 
  theEvaporation = new G4Evaporation;
  theFermiBreakUp = new G4FermiBreakUp;
  theMF = new G4StatMF;
  // evaporation logic
  theHandler = new G4ExcitationHandler;
  theHandler->SetEvaporation(theEvaporation);
  theHandler->SetFermiModel(theFermiBreakUp);
  theHandler->SetMultiFragmentation(theMF);
  theHandler->SetMaxAandZForFermiBreakUp(12, 6);
  theHandler->SetMinEForMultiFrag(3.*MeV);

  // pre-equilibrium stage 
  thePreEquilib = new G4PreCompoundModel(theHandler);
  thePreEquilib->SetMaxEnergy(70*MeV);

  // a no-cascade generator-precompound interaface
  theCascade = new G4GeneratorPrecompoundInterface;
  theCascade->SetDeExcitation(thePreEquilib);

  // QGSP model
  theStringModel = new G4QGSModel<G4QGSParticipants>;
  theHEModel->SetTransport(theCascade);
  theHEModel->SetHighEnergyGenerator(theStringModel);
  theHEModel->SetMinEnergy(6*GeV);
  theHEModel->SetMaxEnergy(100*TeV);
  // Binary cascade for p, n
  theCasc = new G4BinaryCascade;
  theCasc->SetMinEnergy(65*MeV);
  theCasc->SetMaxEnergy(6.1*GeV);
  // fragmentation
  theFragmentation = new G4QGSMFragmentation;
  theStringDecay = new G4ExcitedStringDecay(theFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);
  //
  // Binary Cascade for Pi
  theCascForPi = new G4BinaryCascade;
  theCascForPi->SetMinEnergy(0*MeV);
  theCascForPi->SetMaxEnergy(1.5*GeV);
  // LEP to fill the gap
  theLEPionPlusInelasticModel = new G4LEPionPlusInelastic();
  theLEPionPlusInelasticModel->SetMinEnergy(1.4*GeV);
  theLEPionPlusInelasticModel->SetMaxEnergy(6.1*GeV);
  theLEPionMinusInelasticModel = new G4LEPionMinusInelastic();
  theLEPionMinusInelasticModel->SetMinEnergy(1.4*GeV);
  theLEPionMinusInelasticModel->SetMaxEnergy(6.1*GeV);


  // *******************************************************//
  // *** preparing inelastic reactions for light nuclei *** // 
  // *******************************************************//
  //
  // binary cascade for light nuclei
  // NOTE: Shen XS only up to 10 GeV/n;
  theIonCascade= new G4BinaryLightIonReaction;
  theIonCascade->SetMinEnergy(80*MeV);
  theIonCascade->SetMaxEnergy(40*GeV);
  theTripathiCrossSection = new G4TripathiCrossSection;
  theShenCrossSection = new G4IonsShenCrossSection;
  //
  // deuteron
  theLEDeuteronInelasticModel = new G4LEDeuteronInelastic();
  theLEDeuteronInelasticModel->SetMaxEnergy(100*MeV);
  //
  // triton
  theLETritonInelasticModel = new G4LETritonInelastic();
  theLETritonInelasticModel->SetMaxEnergy(100*MeV);
  //
  // alpha
  theLEAlphaInelasticModel = new G4LEAlphaInelastic();
  theLEAlphaInelasticModel->SetMaxEnergy(100*MeV);
  //
  // Generic Ion and He3
  // NOTE: Shen XS only up to 10 GeV/n;
  theGenIonCascade = new G4BinaryLightIonReaction;
  theGenIonCascade->SetMinEnergy(0*MeV);
  theGenIonCascade->SetMaxEnergy(30*GeV);

  
  // ***************************//
  // *** elastic scattering *** //
  // ***************************//
  //
  theElasticModel = new G4LElastic();
  theElasticProcess.RegisterMe(theElasticModel);


  // *****************************************//
  // *** attaching processes to particles *** //
  // *****************************************//
  //
  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    
    if (particleName == "pi+") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      // NOTE: PreCo crahes for Pi+
      // thePionPlusInelasticProcess.RegisterMe(thePreEquilib);
      thePionPlusInelasticProcess.RegisterMe(theCascForPi);
      thePionPlusInelasticProcess.RegisterMe(theLEPionPlusInelasticModel);
      thePionPlusInelasticProcess.RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(&thePionPlusInelasticProcess);

    } else if (particleName == "pi-") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      // thePionMinusInelasticProcess.RegisterMe(thePreEquilib);
      thePionMinusInelasticProcess.RegisterMe(theCascForPi);
      thePionMinusInelasticProcess.RegisterMe(theLEPionMinusInelasticModel);
      thePionMinusInelasticProcess.RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(&thePionMinusInelasticProcess);
      pmanager->AddRestProcess(&thePiMinusAbsorptionAtRest, ordDefault);

    } else if (particleName == "kaon+") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLEKaonPlusInelasticModel = new G4LEKaonPlusInelastic();
      theLEKaonPlusInelasticModel->SetMaxEnergy(25*GeV);
      theKaonPlusInelasticProcess.RegisterMe(theLEKaonPlusInelasticModel);
      theKaonPlusInelasticProcess.RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(&theKaonPlusInelasticProcess);
      
    } else if (particleName == "kaon0S") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLEKaonZeroSInelasticModel = new G4LEKaonZeroSInelastic();
      theLEKaonZeroSInelasticModel->SetMaxEnergy(25*GeV);
      theKaonZeroSInelasticProcess.RegisterMe(theLEKaonZeroSInelasticModel);
      theKaonZeroSInelasticProcess.RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(&theKaonZeroSInelasticProcess);
      
    } else if (particleName == "kaon0L") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLEKaonZeroLInelasticModel = new G4LEKaonZeroLInelastic();
      theLEKaonZeroLInelasticModel->SetMaxEnergy(25*GeV);
      theKaonZeroLInelasticProcess.RegisterMe(theLEKaonZeroLInelasticModel);
      theKaonZeroLInelasticProcess.RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(&theKaonZeroLInelasticProcess);
      
    } else if (particleName == "kaon-") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLEKaonMinusInelasticModel = new G4LEKaonMinusInelastic();   
      theLEKaonMinusInelasticModel->SetMaxEnergy(25*GeV);
      theKaonMinusInelasticProcess.RegisterMe(theLEKaonMinusInelasticModel);
      theKaonMinusInelasticProcess.RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(&theKaonMinusInelasticProcess);
      pmanager->AddRestProcess(&theKaonMinusAbsorptionAtRest, ordDefault);

    } else if (particleName == "proton") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theProtonInelasticProcess.RegisterMe(thePreEquilib);
      theProtonInelasticProcess.RegisterMe(theCasc);
      theProtonInelasticProcess.RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(&theProtonInelasticProcess);

    } else if (particleName == "anti_proton") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLEAntiProtonInelasticModel = new G4LEAntiProtonInelastic();
      theHEAntiProtonInelasticModel = new G4HEAntiProtonInelastic();
      theLEAntiProtonInelasticModel->SetMaxEnergy(25*GeV);
      theAntiProtonInelasticProcess.RegisterMe(theLEAntiProtonInelasticModel);
      theAntiProtonInelasticProcess.RegisterMe(theHEAntiProtonInelasticModel);
      pmanager->AddDiscreteProcess(&theAntiProtonInelasticProcess);
      pmanager->AddRestProcess(&theAntiProtonAnnihilationAtRest, ordDefault);

    } else if (particleName == "neutron") {
      // elastic scattering
      // LEP
      theNeutronElasticModel1 = new G4LElastic();
      //   theNeutronElasticModel1->SetMinEnergy(19*MeV);
      theNeutronElasticProcess.RegisterMe(theNeutronElasticModel1);
      //   // HP
      //   theNeutronElasticModel2 = new G4NeutronHPElastic();
      //   theNeutronElasticModel2->SetMaxEnergy(19.1*MeV);
      //   theNeutronElasticData = new G4NeutronHPElasticData;
      //   theNeutronElasticProcess.AddDataSet(theNeutronElasticData);
      //   theNeutronElasticProcess.RegisterMe(theNeutronElasticModel2);
      pmanager->AddDiscreteProcess(&theNeutronElasticProcess);
      // inelastic scattering
      //   // HP
      //   theNeutronInelasticModel1 = new G4NeutronHPInelastic();
      //   theNeutronInelasticProcess.RegisterMe(theNeutronInelasticModel1);
      //   theNeutronInelasticData1 = new G4NeutronHPInelasticData;
      //   theNeutronInelasticProcess.AddDataSet(theNeutronInelasticData1);
      // Preco_n + BiC + QGSP
      theNeutronInelasticProcess.RegisterMe(thePreEquilib);
      theNeutronInelasticProcess.RegisterMe(theCasc);
      theNeutronInelasticProcess.RegisterMe(theHEModel);
      pmanager->AddDiscreteProcess(&theNeutronInelasticProcess);
      // capture
      theNeutronCaptureModel1 = new G4LCapture();
      //   theNeutronCaptureModel1->SetMinEnergy(19*MeV);
      theNeutronCaptureProcess.RegisterMe(theNeutronCaptureModel1);
      //   theNeutronCaptureModel2 = new G4NeutronHPCapture;
      //   theNeutronCaptureProcess.RegisterMe(theNeutronCaptureModel2);
      //   theNeutronCaptureData = new G4NeutronHPCaptureData;
      //   theNeutronCaptureProcess.AddDataSet(theNeutronCaptureData);
      pmanager->AddDiscreteProcess(&theNeutronCaptureProcess);
      // fission
      theNeutronFissionModel = new G4LFission();
      theNeutronFissionProcess.RegisterMe(theNeutronFissionModel);
      pmanager->AddDiscreteProcess(&theNeutronFissionProcess);

    } else if (particleName == "anti_neutron") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLEAntiNeutronInelasticModel = new G4LEAntiNeutronInelastic();
      theHEAntiNeutronInelasticModel = new G4HEAntiNeutronInelastic();
      theLEAntiNeutronInelasticModel->SetMaxEnergy(25*GeV);
      theAntiNeutronInelasticProcess.RegisterMe
	(theLEAntiNeutronInelasticModel);
      theAntiNeutronInelasticProcess.RegisterMe
	(theHEAntiNeutronInelasticModel);
      pmanager->AddDiscreteProcess(&theAntiNeutronInelasticProcess);
      pmanager->AddRestProcess(&theAntiNeutronAnnihilationAtRest,ordDefault);

  } else if (particleName == "deuteron") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theDeuteronInelasticProcess = new G4DeuteronInelasticProcess;
      theDeuteronInelasticProcess->AddDataSet(theTripathiCrossSection);
      theDeuteronInelasticProcess->AddDataSet(theShenCrossSection);
      theDeuteronInelasticProcess->RegisterMe(theLEDeuteronInelasticModel);
      theDeuteronInelasticProcess->RegisterMe(theIonCascade);
      pmanager->AddDiscreteProcess(theDeuteronInelasticProcess);

    } else if (particleName == "triton") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theTritonInelasticProcess = new G4TritonInelasticProcess;
      theTritonInelasticProcess->AddDataSet(theTripathiCrossSection);
      theTritonInelasticProcess->AddDataSet(theShenCrossSection);
      theTritonInelasticProcess->RegisterMe(theLETritonInelasticModel);
      theTritonInelasticProcess->RegisterMe(theIonCascade);
      pmanager->AddDiscreteProcess(theTritonInelasticProcess);

    } else if (particleName == "alpha") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theAlphaInelasticProcess = new G4AlphaInelasticProcess;
      theAlphaInelasticProcess->AddDataSet(theTripathiCrossSection);
      theAlphaInelasticProcess->AddDataSet(theShenCrossSection);
      theAlphaInelasticProcess->RegisterMe(theLEAlphaInelasticModel);
      theAlphaInelasticProcess->RegisterMe(theIonCascade);
      pmanager->AddDiscreteProcess(theAlphaInelasticProcess);
      
    } else if (particleName == "He3") {
      // NOTE elastic scattering does not stick to He3!
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theHe3InelasticProcess = new G4HadronInelasticProcess
	("He3Inelastic", G4He3::He3());
      theHe3InelasticProcess->AddDataSet(theTripathiCrossSection);
      theHe3InelasticProcess->AddDataSet(theShenCrossSection);
      theHe3InelasticProcess->RegisterMe(theGenIonCascade);
      pmanager->AddDiscreteProcess(theHe3InelasticProcess);

    } else if (particleName == "GenericIon") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theGenericIonInelasticProcess = new G4HadronInelasticProcess
	("IonInelastic", G4GenericIon::GenericIon());
      theGenericIonInelasticProcess->AddDataSet(theTripathiCrossSection);
      theGenericIonInelasticProcess->AddDataSet(theShenCrossSection);
      theGenericIonInelasticProcess->RegisterMe(theGenIonCascade);
      pmanager->AddDiscreteProcess(theGenericIonInelasticProcess);
      
    } else if (particleName == "lambda") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLELambdaInelasticModel = new G4LELambdaInelastic();
      theHELambdaInelasticModel = new G4HELambdaInelastic();
      theLELambdaInelasticModel->SetMaxEnergy(25*GeV);
      theLambdaInelasticProcess.RegisterMe(theLELambdaInelasticModel);
      theLambdaInelasticProcess.RegisterMe(theHELambdaInelasticModel);
      pmanager->AddDiscreteProcess(&theLambdaInelasticProcess);

    } else if (particleName == "anti_lambda") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLEAntiLambdaInelasticModel = new G4LEAntiLambdaInelastic();
      theHEAntiLambdaInelasticModel = new G4HEAntiLambdaInelastic();
      theLEAntiLambdaInelasticModel->SetMaxEnergy(25*GeV);
      theAntiLambdaInelasticProcess.RegisterMe(theLEAntiLambdaInelasticModel);
      theAntiLambdaInelasticProcess.RegisterMe(theHEAntiLambdaInelasticModel);
      pmanager->AddDiscreteProcess(&theAntiLambdaInelasticProcess);

    } else if (particleName == "omega-") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLEOmegaMinusInelasticModel = new G4LEOmegaMinusInelastic();
      theHEOmegaMinusInelasticModel = new G4HEOmegaMinusInelastic();
      theLEOmegaMinusInelasticModel->SetMaxEnergy(25*GeV);
      theOmegaMinusInelasticProcess.RegisterMe(theLEOmegaMinusInelasticModel);
      theOmegaMinusInelasticProcess.RegisterMe(theHEOmegaMinusInelasticModel);
      pmanager->AddDiscreteProcess(&theOmegaMinusInelasticProcess);

    } else if (particleName == "anti_omega-") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLEAntiOmegaMinusInelasticModel = new G4LEAntiOmegaMinusInelastic(); 
      theHEAntiOmegaMinusInelasticModel = new G4HEAntiOmegaMinusInelastic(); 
      theLEAntiOmegaMinusInelasticModel->SetMaxEnergy(25*GeV);
      theAntiOmegaMinusInelasticProcess.RegisterMe
	(theLEAntiOmegaMinusInelasticModel);
      theAntiOmegaMinusInelasticProcess.RegisterMe
	(theHEAntiOmegaMinusInelasticModel);
      pmanager->AddDiscreteProcess(&theAntiOmegaMinusInelasticProcess);

    } else if (particleName == "sigma-") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLESigmaMinusInelasticModel = new G4LESigmaMinusInelastic();
      theHESigmaMinusInelasticModel = new G4HESigmaMinusInelastic();
      theLESigmaMinusInelasticModel->SetMaxEnergy(25*GeV);
      theSigmaMinusInelasticProcess.RegisterMe(theLESigmaMinusInelasticModel);
      theSigmaMinusInelasticProcess.RegisterMe(theHESigmaMinusInelasticModel);
      pmanager->AddDiscreteProcess(&theSigmaMinusInelasticProcess);

    } else if (particleName == "anti_sigma-") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLEAntiSigmaMinusInelasticModel = new G4LEAntiSigmaMinusInelastic();
      theHEAntiSigmaMinusInelasticModel = new G4HEAntiSigmaMinusInelastic();
      theLEAntiSigmaMinusInelasticModel->SetMaxEnergy(25*GeV);
      theAntiSigmaMinusInelasticProcess.RegisterMe
	(theLEAntiSigmaMinusInelasticModel);
      theAntiSigmaMinusInelasticProcess.RegisterMe
	(theHEAntiSigmaMinusInelasticModel);
      pmanager->AddDiscreteProcess(&theAntiSigmaMinusInelasticProcess);

    } else if (particleName == "sigma+") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLESigmaPlusInelasticModel = new G4LESigmaPlusInelastic();
      theHESigmaPlusInelasticModel = new G4HESigmaPlusInelastic();      
      theLESigmaPlusInelasticModel->SetMaxEnergy(25*GeV);
      theSigmaPlusInelasticProcess.RegisterMe(theLESigmaPlusInelasticModel);
      theSigmaPlusInelasticProcess.RegisterMe(theHESigmaPlusInelasticModel);
      pmanager->AddDiscreteProcess(&theSigmaPlusInelasticProcess);

    } else if (particleName == "anti_sigma+") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLEAntiSigmaPlusInelasticModel = new G4LEAntiSigmaPlusInelastic();
      theHEAntiSigmaPlusInelasticModel = new G4HEAntiSigmaPlusInelastic();
      theLEAntiSigmaPlusInelasticModel->SetMaxEnergy(25*GeV);
      theAntiSigmaPlusInelasticProcess.RegisterMe
	(theLEAntiSigmaPlusInelasticModel);
      theAntiSigmaPlusInelasticProcess.RegisterMe
	(theHEAntiSigmaPlusInelasticModel);
      pmanager->AddDiscreteProcess(&theAntiSigmaPlusInelasticProcess);

    } else if (particleName == "xi0") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLEXiZeroInelasticModel = new G4LEXiZeroInelastic();
      theHEXiZeroInelasticModel = new G4HEXiZeroInelastic();
      theLEXiZeroInelasticModel->SetMaxEnergy(25*GeV);
      theXiZeroInelasticProcess.RegisterMe(theLEXiZeroInelasticModel);
      theXiZeroInelasticProcess.RegisterMe(theHEXiZeroInelasticModel);
      pmanager->AddDiscreteProcess(&theXiZeroInelasticProcess);

    } else if (particleName == "anti_xi0") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLEAntiXiZeroInelasticModel = new G4LEAntiXiZeroInelastic();
      theHEAntiXiZeroInelasticModel = new G4HEAntiXiZeroInelastic();
      theLEAntiXiZeroInelasticModel->SetMaxEnergy(25*GeV);
      theAntiXiZeroInelasticProcess.RegisterMe(theLEAntiXiZeroInelasticModel);
      theAntiXiZeroInelasticProcess.RegisterMe(theHEAntiXiZeroInelasticModel);
      pmanager->AddDiscreteProcess(&theAntiXiZeroInelasticProcess);

    } else if (particleName == "xi-") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLEXiMinusInelasticModel = new G4LEXiMinusInelastic();
      theHEXiMinusInelasticModel = new G4HEXiMinusInelastic();
      theLEXiMinusInelasticModel->SetMaxEnergy(25*GeV);
      theXiMinusInelasticProcess.RegisterMe(theLEXiMinusInelasticModel);
      theXiMinusInelasticProcess.RegisterMe(theHEXiMinusInelasticModel);
      pmanager->AddDiscreteProcess(&theXiMinusInelasticProcess);

    } else if (particleName == "anti_xi-") {
      pmanager->AddDiscreteProcess(&theElasticProcess);
      theLEAntiXiMinusInelasticModel = new G4LEAntiXiMinusInelastic();
      theHEAntiXiMinusInelasticModel = new G4HEAntiXiMinusInelastic();
      theLEAntiXiMinusInelasticModel->SetMaxEnergy(25*GeV);
      theAntiXiMinusInelasticProcess.RegisterMe
	(theLEAntiXiMinusInelasticModel);
      theAntiXiMinusInelasticProcess.RegisterMe
	(theHEAntiXiMinusInelasticModel);
      pmanager->AddDiscreteProcess(&theAntiXiMinusInelasticProcess);
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

