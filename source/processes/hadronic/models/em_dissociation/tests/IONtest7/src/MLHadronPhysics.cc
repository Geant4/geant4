////////////////////////////////////////////////////////////////////////////////
//
#include "MLHadronPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>
////////////////////////////////////////////////////////////////////////////////
//
MLHadronPhysics::MLHadronPhysics (const G4String& name)
  : G4VPhysicsConstructor(name), mode(name)
{}
////////////////////////////////////////////////////////////////////////////////
//
MLHadronPhysics::~MLHadronPhysics ()
{}
////////////////////////////////////////////////////////////////////////////////
//
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

// Nuclei
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
////////////////////////////////////////////////////////////////////////////////
//
void MLHadronPhysics::ConstructParticle ()
{
  //  Construct all mesons
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  //  Construct all barions
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  

}
////////////////////////////////////////////////////////////////////////////////
//
#include "G4ProcessManager.hh"

// -- geneator models
#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4Evaporation.hh"
#include "G4CompetitiveFission.hh"
#include "G4FermiBreakUp.hh"
#include "G4StatMF.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4LEProtonInelastic.hh"
#include "G4StringModel.hh"
#include "G4PreCompoundModel.hh"
#include "G4FTFModel.hh"
#include "G4QGSModel.hh"
#include "G4QGSParticipants.hh"
#include "G4LundStringFragmentation.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4StringChipsParticleLevelInterface.hh"

//HPNeutron

#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"

#include "G4Mars5GeV.hh"
#include "G4BinaryCascade.hh"
#include "G4CascadeInterface.hh"

#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"

void MLHadronPhysics::ConstructProcess()
{
  G4ProcessManager * pManager = 0;
  
  // first to decide what prcesses shall be included according to 
  // the mode string passed on
  //
  //
  G4bool MARS5GEV = false;
  G4bool HPNEUTRON = true;
  G4bool PRECOMP = false; // was true before 22/09/03
  G4bool BC = true ;    // from 22/09/03 the default cascade code is changed to BC
  G4bool CC = false ;

  // default mode is LE+PC+LN (Hadron)
   if ( mode == "Hadron-n" ) {
    HPNEUTRON = false ;
  }else if ( mode == "Hadron-pc-n" ) {
    HPNEUTRON = false ;
    PRECOMP = false;
  } else if ( mode == "Mars5GeV" ) {
    MARS5GEV = true;
    HPNEUTRON = true ;
    PRECOMP = false;
  } else if ( mode == "Mars5GeV+pc" ) {
    MARS5GEV = true;
    PRECOMP = true; 
    HPNEUTRON = true ;
  } else if ( mode == "Mars5GeV-ln" ) {
    MARS5GEV = true;
    PRECOMP = false;
    HPNEUTRON = false;
  } else if ( mode == "Binary") {
    BC = true ;
    CC = false;
    PRECOMP = false; 
    HPNEUTRON = true ;   
  } else if ( mode == "Classic") {
    CC = true ;
    BC = false ;
    PRECOMP = false; 
    HPNEUTRON = true ;
  }
  // this will be the model class for high energies
  G4TheoFSGenerator * theTheoModel = new G4TheoFSGenerator;
   
  // all models for treatment of thermal nucleus 
  G4Evaporation * theEvaporation = new G4Evaporation;
  G4FermiBreakUp * theFermiBreakUp = new G4FermiBreakUp;
  G4StatMF * theMF = new G4StatMF;

  // Evaporation logic
  G4ExcitationHandler * theHandler = new G4ExcitationHandler;
  theHandler->SetEvaporation(theEvaporation);
  theHandler->SetFermiModel(theFermiBreakUp);
  theHandler->SetMultiFragmentation(theMF);
  theHandler->SetMaxAandZForFermiBreakUp(12, 6);
  theHandler->SetMinEForMultiFrag(5*MeV);
	
  // Pre equilibrium stage 
  G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel(theHandler);
    
  // a no-cascade generator-precompound interaface
  G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface;
  theCascade->SetDeExcitation(thePreEquilib);  
	
  // a cascade interface interface to chiral invriant phase space decay
  // at particle level.
  //  G4StringChipsParticleLevelInterface * theCascade = new G4StringChipsParticleLevelInterface();

  // here come the high energy parts
  // the string model; still not quite according to design - Explicite use of the forseen interfaces     
  G4VPartonStringModel * theStringModel;
  theStringModel = new G4QGSModel<G4QGSParticipants>;
  theTheoModel->SetTransport(theCascade);
  theTheoModel->SetHighEnergyGenerator(theStringModel);
  theTheoModel->SetMinEnergy(9*GeV);  // 15 GeV may be the right limit
  theTheoModel->SetMaxEnergy(100*TeV);
  
  //  G4VLongitudinalStringDecay * theFragmentation = new G4LundStringFragmentation;
  G4VLongitudinalStringDecay * theFragmentation = new G4QGSMFragmentation;
  G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay(theFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);
  
 // G4VPartonStringModel * theStringModel = new G4FTFModel;
 // theTheoModel->SetTransport(theCascade);
 // theTheoModel->SetHighEnergyGenerator(theStringModel);
 // theTheoModel->SetMinEnergy(15*GeV);
 // theTheoModel->SetMaxEnergy(100*TeV);
    
//  G4VLongitudinalStringDecay * theFragmentation = new G4QGSMFragmentation;
//  G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay(theFragmentation);
//  theStringModel->SetFragmentationModel(theStringDecay);

  // PreCompound
  G4PreCompoundModel * thePreCompoundModel = new G4PreCompoundModel(theHandler);
  thePreCompoundModel->SetMaxEnergy(120.*MeV);
  
  // Elastic Process
  theElasticModel = new G4LElastic();
  theElasticProcess.RegisterMe(theElasticModel);

  // Mars5GeV model
  G4Mars5GeV* theMars5GeVModel = new G4Mars5GeV();
  theMars5GeVModel->SetMaxEnergy(5.*GeV);

  // ---------------------------------------------------------------------------

  // PionPlus
  pManager = G4PionPlus::PionPlus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEPionPlusModel = new G4LEPionPlusInelastic();
  
  if (MARS5GEV) {
    thePionPlusInelastic.RegisterMe(theMars5GeVModel);
    theLEPionPlusModel->SetMaxEnergy(1.*MeV);    
  }
  thePionPlusInelastic.RegisterMe(theLEPionPlusModel);
  
  theHEPionPlusModel = new G4HEPionPlusInelastic();
  thePionPlusInelastic.RegisterMe(theHEPionPlusModel);
  pManager->AddDiscreteProcess(&thePionPlusInelastic);
  
  // PionMinus
  pManager = G4PionMinus::PionMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEPionMinusModel = new G4LEPionMinusInelastic();
  if (MARS5GEV) {
    thePionMinusInelastic.RegisterMe(theMars5GeVModel);
    theLEPionMinusModel->SetMaxEnergy(1.*MeV);  
  }
  thePionMinusInelastic.RegisterMe(theLEPionMinusModel);
  
  theHEPionMinusModel = new G4HEPionMinusInelastic();
  thePionMinusInelastic.RegisterMe(theHEPionMinusModel);
  pManager->AddDiscreteProcess(&thePionMinusInelastic);

  pManager->AddRestProcess(&thePionMinusAbsorption, ordDefault);

  // KaonPlus
  pManager = G4KaonPlus::KaonPlus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEKaonPlusModel = new G4LEKaonPlusInelastic();
  if (MARS5GEV) {
    theKaonPlusInelastic.RegisterMe(theMars5GeVModel);
    theLEKaonPlusModel->SetMaxEnergy(1.*MeV); 
  } 
  theKaonPlusInelastic.RegisterMe(theLEKaonPlusModel);
  
  theHEKaonPlusModel = new G4HEKaonPlusInelastic();
  theKaonPlusInelastic.RegisterMe(theHEKaonPlusModel);
  pManager->AddDiscreteProcess(&theKaonPlusInelastic);

  // KaonMinus
  pManager = G4KaonMinus::KaonMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEKaonMinusModel = new G4LEKaonMinusInelastic();
  if (MARS5GEV) {
    theKaonMinusInelastic.RegisterMe(theMars5GeVModel);
    theLEKaonMinusModel->SetMaxEnergy(1.*MeV);
  }
  theKaonMinusInelastic.RegisterMe(theLEKaonMinusModel);
  
  theHEKaonMinusModel = new G4HEKaonMinusInelastic();
  theKaonMinusInelastic.RegisterMe(theHEKaonMinusModel);
  pManager->AddDiscreteProcess(&theKaonMinusInelastic);

  // KaonZeroL
  pManager = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEKaonZeroLModel = new G4LEKaonZeroLInelastic();
  theHEKaonZeroLModel = new G4HEKaonZeroInelastic();
  theKaonZeroLInelastic.RegisterMe(theLEKaonZeroLModel);
  theKaonZeroLInelastic.RegisterMe(theHEKaonZeroLModel);
  pManager->AddDiscreteProcess(&theKaonZeroLInelastic);
 
  // KaonZeroS
  pManager = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEKaonZeroSModel = new G4LEKaonZeroSInelastic();
  theHEKaonZeroSModel = new G4HEKaonZeroInelastic();
  theKaonZeroSInelastic.RegisterMe(theLEKaonZeroSModel);
  theKaonZeroSInelastic.RegisterMe(theHEKaonZeroSModel);
  pManager->AddDiscreteProcess(&theKaonZeroSInelastic);

  // Proton
  pManager = G4Proton::Proton()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);
  theLEProtonModel = new G4LEProtonInelastic();

  if (MARS5GEV) {
    G4Mars5GeV* theMars5GeVModelforProton = new G4Mars5GeV();
    theMars5GeVModelforProton->SetMaxEnergy(5.*GeV);
    if (PRECOMP) {
      theMars5GeVModelforProton->SetMinEnergy(120.*MeV);
      theProtonInelastic.RegisterMe(thePreCompoundModel);
    } else {
      theMars5GeVModelforProton->SetMinEnergy(1.*MeV);
    }
    theLEProtonModel->SetMaxEnergy(1*MeV);
    theProtonInelastic.RegisterMe(theLEProtonModel);
    theProtonInelastic.RegisterMe(theMars5GeVModelforProton);
  } else if (BC) {
    // Binary Cascade
    G4BinaryCascade * theBC = new G4BinaryCascade;
    theBC->SetMaxEnergy(10.*GeV);
    theProtonInelastic.RegisterMe(theBC);
    //    theLEProtonModel->SetMinEnergy(9.5*GeV);
    //   theLEProtonModel->SetMaxEnergy(15.5*GeV);
    theProtonInelastic.RegisterMe(theTheoModel);
  } else if (CC) {
    // Classis Cascade
    G4CascadeInterface * theCC = new G4CascadeInterface;
    theCC->SetMaxEnergy(10.*GeV);
    theProtonInelastic.RegisterMe(theCC);
    //    theLEProtonModel->SetMinEnergy(9.5*GeV);
    // theLEProtonModel->SetMaxEnergy(15.5*GeV);
    theProtonInelastic.RegisterMe(theTheoModel); 
  } else {
    if (PRECOMP) {
      theLEProtonModel->SetMinEnergy(100.*MeV);
      G4PreCompoundModel* thePreCompoundModel = new G4PreCompoundModel(theHandler);
      thePreCompoundModel->SetMaxEnergy(110.*MeV);
      theProtonInelastic.RegisterMe(thePreCompoundModel);
    }
    theProtonInelastic.RegisterMe(theLEProtonModel);
    theHEProtonModel = new G4HEProtonInelastic();
    theProtonInelastic.RegisterMe(theHEProtonModel);
  }
  // now the cross-sections.
  //  G4CrossSectionDataStore * theStoreProton =
  //  ((G4HadronInelasticProcess*)&theProtonInelastic)->GetCrossSectionDataStore();
  G4ProtonInelasticCrossSection * theProtonData = new G4ProtonInelasticCrossSection;
  //theStoreProton->AddDataSet(theProtonData);
  theProtonInelastic.AddDataSet(theProtonData);
  pManager->AddDiscreteProcess(&theProtonInelastic);

  // anti-Proton
  pManager = G4AntiProton::AntiProton()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);
  
  theLEAntiProtonModel = new G4LEAntiProtonInelastic();
  theHEAntiProtonModel = new G4HEAntiProtonInelastic();
  theAntiProtonInelastic.RegisterMe(theLEAntiProtonModel);
  theAntiProtonInelastic.RegisterMe(theHEAntiProtonModel);
  pManager->AddDiscreteProcess(&theAntiProtonInelastic);

  pManager->AddRestProcess(&theAntiProtonAnnihilation);

  // Neutron
  pManager = G4Neutron::Neutron()->GetProcessManager();
  // add process

  // elastic scattering
  if (HPNEUTRON) {
    G4HadronElasticProcess* theNeutronElasticProcess = 
      new G4HadronElasticProcess();
    G4LElastic* theElasticModel1 = new G4LElastic;
    G4NeutronHPElastic * theElasticNeutron = new G4NeutronHPElastic;
    theNeutronElasticProcess->RegisterMe(theElasticModel1);
    theElasticModel1->SetMinEnergy(19.*MeV);
    theNeutronElasticProcess->RegisterMe(theElasticNeutron);
    theElasticNeutron->SetMaxEnergy(20.*MeV);
    //   G4CrossSectionDataStore * theStore = 
    //   ((G4HadronElasticProcess*)theNeutronElasticProcess)
    //  ->GetCrossSectionDataStore();
    G4NeutronHPElasticData * theNeutronData = new G4NeutronHPElasticData;
    //    theStore->AddDataSet(theNeutronData);
    theNeutronElasticProcess->AddDataSet(theNeutronData);
    pManager->AddDiscreteProcess(theNeutronElasticProcess);
  } else {
    pManager->AddDiscreteProcess(&theElasticProcess);
  }

  // inelastic 
  G4Mars5GeV* theMars5GeVModelforNeutron = new G4Mars5GeV();   
  theMars5GeVModelforNeutron->SetMaxEnergy(5.*GeV);
  if (MARS5GEV) {
    theNeutronInelastic.RegisterMe(theMars5GeVModelforNeutron);
  }
  G4NeutronHPInelastic * theHPNeutronInelasticModel =
    new G4NeutronHPInelastic;
  theHPNeutronInelasticModel->SetMaxEnergy(20.*MeV);
  if (HPNEUTRON) {
    theNeutronInelastic.RegisterMe(theHPNeutronInelasticModel);
    //   G4CrossSectionDataStore * theStore1 = 
    //  ((G4HadronInelasticProcess*) &theNeutronInelastic)
    //  ->GetCrossSectionDataStore();
    G4NeutronHPInelasticData * theNeutronData1 = new G4NeutronHPInelasticData;
    theNeutronInelastic.AddDataSet(theNeutronData1);
  }
  G4PreCompoundModel* neutronPreCompoundModel = new G4PreCompoundModel(theHandler);
  neutronPreCompoundModel->SetMaxEnergy(120.*MeV);
  if ( PRECOMP) {
    theNeutronInelastic.RegisterMe(neutronPreCompoundModel);
  }
  
  G4LENeutronInelastic* neutronLEInelasticModel = new G4LENeutronInelastic;
  //theNeutronInelastic.RegisterMe(neutronLEInelasticModel);

  //
  if (MARS5GEV ) {
    neutronLEInelasticModel->SetMinEnergy(5.*GeV);
    if (PRECOMP) { 
      theMars5GeVModelforNeutron->SetMinEnergy(100.*MeV);
      if (HPNEUTRON) neutronPreCompoundModel->SetMinEnergy(19.*MeV);
    } else {
      if (HPNEUTRON) theMars5GeVModelforNeutron->SetMinEnergy(19.*MeV);
    }
    theNeutronInelastic.RegisterMe(neutronLEInelasticModel);
  } else if (BC) {    
    G4BinaryCascade * neutronBC = new G4BinaryCascade;
    if (HPNEUTRON) neutronBC->SetMinEnergy(19.*MeV);
    neutronBC->SetMaxEnergy(10.*GeV);
    theNeutronInelastic.RegisterMe(neutronBC);
    //    neutronLEInelasticModel->SetMinEnergy(9.5*GeV);
    //neutronLEInelasticModel->SetMaxEnergy(15.5*GeV);
    theNeutronInelastic.RegisterMe(theTheoModel);  
  } else if (CC) {
    G4CascadeInterface * neutronCC = new G4CascadeInterface;
    if (HPNEUTRON) neutronCC->SetMinEnergy(19.*MeV);
    neutronCC->SetMaxEnergy(10.*GeV);
    theNeutronInelastic.RegisterMe(neutronCC);
    //neutronLEInelasticModel->SetMinEnergy(9.5*GeV);
    //neutronLEInelasticModel->SetMaxEnergy(15.5*GeV);
    theNeutronInelastic.RegisterMe(theTheoModel);  
  } else {
    if (PRECOMP) { 
      neutronLEInelasticModel->SetMinEnergy(100.*MeV);
      if (HPNEUTRON) neutronPreCompoundModel->SetMinEnergy(19.*MeV);
    } else {
      if (HPNEUTRON) { 
	neutronLEInelasticModel->SetMinEnergy(19.*MeV);
      } 
    }  
    theNeutronInelastic.RegisterMe(neutronLEInelasticModel);
    theHENeutronModel = new G4HENeutronInelastic();
    theNeutronInelastic.RegisterMe(theHENeutronModel);
  }  
  
  // now the cross-sections.
  //  G4CrossSectionDataStore * theStoreNeutron =
  //  ((G4HadronInelasticProcess*)&theNeutronInelastic)->GetCrossSectionDataStore();
  G4NeutronInelasticCrossSection * theNeutronData = new G4NeutronInelasticCrossSection;
  //theStoreNeutron->AddDataSet(theNeutronData);
  theNeutronInelastic.AddDataSet(theNeutronData);
  pManager->AddDiscreteProcess(&theNeutronInelastic);
  
  // fission
  G4HadronFissionProcess* theFissionProcess =
    new G4HadronFissionProcess;
  G4LFission* theFissionModel = new G4LFission;
  theFissionProcess->RegisterMe(theFissionModel);
  pManager->AddDiscreteProcess(theFissionProcess);
  
  //capture  
  G4HadronCaptureProcess* theCaptureProcess =
    new G4HadronCaptureProcess;
  G4LCapture* theCaptureModel = new G4LCapture;
  theCaptureProcess->RegisterMe(theCaptureModel);
  if (HPNEUTRON) {
    theCaptureModel->SetMinEnergy(19.*MeV);
    G4NeutronHPCapture * theLENeutronCaptureModel = new G4NeutronHPCapture;
    theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
    //    G4CrossSectionDataStore * theStore3 = 
    //  ((G4HadronCaptureProcess*)theCaptureProcess)->
    //  GetCrossSectionDataStore();
    G4NeutronHPCaptureData * theNeutronData3 = new G4NeutronHPCaptureData;
    //theStore3->AddDataSet(theNeutronData3);
    theCaptureProcess->AddDataSet(theNeutronData3);
  }
  pManager->AddDiscreteProcess(theCaptureProcess);
  
  // AntiNeutron
  pManager = G4AntiNeutron::AntiNeutron()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);
  
  theLEAntiNeutronModel = new G4LEAntiNeutronInelastic();
  theHEAntiNeutronModel = new G4HEAntiNeutronInelastic();
  theAntiNeutronInelastic.RegisterMe(theLEAntiNeutronModel);
  theAntiNeutronInelastic.RegisterMe(theHEAntiNeutronModel);
  pManager->AddDiscreteProcess(&theAntiNeutronInelastic);
    
  pManager->AddRestProcess(&theAntiNeutronAnnihilation);

  // Lambda
  pManager = G4Lambda::Lambda()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLELambdaModel = new G4LELambdaInelastic();
  theHELambdaModel = new G4HELambdaInelastic();
  theLambdaInelastic.RegisterMe(theLELambdaModel);
  theLambdaInelastic.RegisterMe(theHELambdaModel);
  pManager->AddDiscreteProcess(&theLambdaInelastic);
  
  // AntiLambda
  pManager = G4AntiLambda::AntiLambda()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);
  
  theLEAntiLambdaModel = new G4LEAntiLambdaInelastic();
  theHEAntiLambdaModel = new G4HEAntiLambdaInelastic();
  theAntiLambdaInelastic.RegisterMe(theLEAntiLambdaModel);
  theAntiLambdaInelastic.RegisterMe(theHEAntiLambdaModel);
  pManager->AddDiscreteProcess(&theAntiLambdaInelastic);
    
  // SigmaMinus
  pManager = G4SigmaMinus::SigmaMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLESigmaMinusModel = new G4LESigmaMinusInelastic();
  theHESigmaMinusModel = new G4HESigmaMinusInelastic();
  theSigmaMinusInelastic.RegisterMe(theLESigmaMinusModel);
  theSigmaMinusInelastic.RegisterMe(theHESigmaMinusModel);
  pManager->AddDiscreteProcess(&theSigmaMinusInelastic);

  // anti-SigmaMinus
  pManager = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiSigmaMinusModel = new G4LEAntiSigmaMinusInelastic();
  theHEAntiSigmaMinusModel = new G4HEAntiSigmaMinusInelastic();
  theAntiSigmaMinusInelastic.RegisterMe(theLEAntiSigmaMinusModel);
  theAntiSigmaMinusInelastic.RegisterMe(theHEAntiSigmaMinusModel);
  pManager->AddDiscreteProcess(&theAntiSigmaMinusInelastic);

  // SigmaPlus
  pManager = G4SigmaPlus::SigmaPlus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLESigmaPlusModel = new G4LESigmaPlusInelastic();
  theHESigmaPlusModel = new G4HESigmaPlusInelastic();
  theSigmaPlusInelastic.RegisterMe(theLESigmaPlusModel);
  theSigmaPlusInelastic.RegisterMe(theHESigmaPlusModel);
  pManager->AddDiscreteProcess(&theSigmaPlusInelastic);

  // anti-SigmaPlus
  pManager = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiSigmaPlusModel = new G4LEAntiSigmaPlusInelastic();
  theHEAntiSigmaPlusModel = new G4HEAntiSigmaPlusInelastic();
  theAntiSigmaPlusInelastic.RegisterMe(theLEAntiSigmaPlusModel);
  theAntiSigmaPlusInelastic.RegisterMe(theHEAntiSigmaPlusModel);
  pManager->AddDiscreteProcess(&theAntiSigmaPlusInelastic);

  // XiMinus
  pManager = G4XiMinus::XiMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);
  
  theLEXiMinusModel = new G4LEXiMinusInelastic();
  theHEXiMinusModel = new G4HEXiMinusInelastic();
  theXiMinusInelastic.RegisterMe(theLEXiMinusModel);
  theXiMinusInelastic.RegisterMe(theHEXiMinusModel);
  pManager->AddDiscreteProcess(&theXiMinusInelastic);

  // anti-XiMinus
  pManager = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiXiMinusModel = new G4LEAntiXiMinusInelastic();
  theHEAntiXiMinusModel = new G4HEAntiXiMinusInelastic();
  theAntiXiMinusInelastic.RegisterMe(theLEAntiXiMinusModel);
  theAntiXiMinusInelastic.RegisterMe(theHEAntiXiMinusModel);
  pManager->AddDiscreteProcess(&theAntiXiMinusInelastic);

  // XiZero
  pManager = G4XiZero::XiZero()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEXiZeroModel = new G4LEXiZeroInelastic();
  theHEXiZeroModel = new G4HEXiZeroInelastic();
  theXiZeroInelastic.RegisterMe(theLEXiZeroModel);
  theXiZeroInelastic.RegisterMe(theHEXiZeroModel);
  pManager->AddDiscreteProcess(&theXiZeroInelastic);

  // anti-XiZero
  pManager = G4AntiXiZero::AntiXiZero()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiXiZeroModel = new G4LEAntiXiZeroInelastic();
  theHEAntiXiZeroModel = new G4HEAntiXiZeroInelastic();
  theAntiXiZeroInelastic.RegisterMe(theLEAntiXiZeroModel);
  theAntiXiZeroInelastic.RegisterMe(theHEAntiXiZeroModel);
  pManager->AddDiscreteProcess(&theAntiXiZeroInelastic);

  // OmegaMinus
  pManager = G4OmegaMinus::OmegaMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEOmegaMinusModel = new G4LEOmegaMinusInelastic();
  theHEOmegaMinusModel = new G4HEOmegaMinusInelastic();
  theOmegaMinusInelastic.RegisterMe(theLEOmegaMinusModel);
  theOmegaMinusInelastic.RegisterMe(theHEOmegaMinusModel);
  pManager->AddDiscreteProcess(&theOmegaMinusInelastic);

  // anti-OmegaMinus
  pManager = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();
  // add process
  pManager->AddDiscreteProcess(&theElasticProcess);

  theLEAntiOmegaMinusModel = new G4LEAntiOmegaMinusInelastic();
  theHEAntiOmegaMinusModel = new G4HEAntiOmegaMinusInelastic();
  theAntiOmegaMinusInelastic.RegisterMe(theLEAntiOmegaMinusModel);
  theAntiOmegaMinusInelastic.RegisterMe(theHEAntiOmegaMinusModel);
  pManager->AddDiscreteProcess(&theAntiOmegaMinusInelastic);

}
////////////////////////////////////////////////////////////////////////////////
