#include "G4DosimPhysics.hh"
#include "Hall.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "Model.hh"
#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4LElastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LEAlphaInelastic.hh"
#include "G4LEDeuteronInelastic.hh"

#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronHPFission.hh"
#include "G4NeutronHPFissionData.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4CrossSectionDataStore.hh"

#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4PreCompoundNeutron.hh"
#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"

#include "G4Mars5GeV.hh"
#include "Mars01EminCut.hh"

#include "G4hIonisation.hh"

#include "G4PionMinusNuclearReaction.hh"

//#include "stdafx.h"
PhysicsList::~PhysicsList()
{
  //Mahame procesite
  if(pLowEnergyModel) delete pLowEnergyModel;
  if(pPrecompoundModel) delete pPrecompoundModel;
  if(pMars5GeVModel) delete pMars5GeVModel;
  if(pChiralInvModel) delete pChiralInvModel;
  if(pProtonLowEnergyModel) delete pProtonLowEnergyModel;
  if(pProtonPrecompoundModel) delete pProtonPrecompoundModel;
  if(pProtonMars5GeVModel) delete pProtonMars5GeVModel;
  if(pCutProcess) delete pCutProcess;
  if(pProtonCutProcess) delete pProtonCutProcess;
  if(pProtonChiralInvModel) delete pProtonChiralInvModel;
  if(pAlphaLowEnergy) delete pAlphaLowEnergy;
  if(pAlphaPrecompound) delete pAlphaPrecompound;
  if(pAlphaMars) delete pAlphaMars;
  if(pAlphaCut) delete pAlphaCut;
  if(pAlphaChiral) delete pAlphaChiral;
  if(pDeuteronLowEnergy) delete pDeuteronLowEnergy;
  if(pDeuteronPrecompound) delete pDeuteronPrecompound;
  if(pDeuteronMars) delete pDeuteronMars;
  if(pDeuteronCut) delete pDeuteronCut;
  if(pDeuteronChiral) delete pDeuteronChiral;
  delete pMessenger;
}
void PhysicsList::ConstructParticle()
{
  if(G4Proton::ProtonDefinition());
  if(G4Neutron::NeutronDefinition());
  G4Alpha::AlphaDefinition();
  G4Deuteron::DeuteronDefinition();
  G4He3::He3Definition();
}
void PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructParticleProcess();
  ConstructGlobal();
}

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "MyLDP.hh"
#include "G4VEvaporationChannel.hh"
#include "EvapChannel.hh"
#include "G4Evaporation.hh"
#include "ProtChannel.hh"
#include "NeutChannel.hh"
#include "DeutChannel.hh"
#include "TritChannel.hh"
#include "He3Channel.hh"
#include "AlphaChannel.hh"
#include "ProtonEvapProb.hh"
#include "NeutronEvapProb.hh"
#include "AlphaEvapProb.hh"
#include "TritonEvapProb.hh"
#include "DeuteronEvapProb.hh"
#include "He3EvapProb.hh"
#include "G4CompetitiveFission.hh"
#include "G4PhotonEvaporation.hh"
#include "FissionLDP.hh"
#include "MyNewClass.hh"

void PhysicsList::ConstructParticleProcess()
{
  G4HadronElasticProcess* pElasticProcess = new G4HadronElasticProcess;
  G4LElastic* pElasticModel = new G4LElastic();
  pElasticProcess->RegisterMe(pElasticModel);
  G4HadronElasticProcess* pElasticProcess1 = new G4HadronElasticProcess;

  theParticleIterator->reset();
  while((*theParticleIterator)()){
    G4ParticleDefinition* pParticle = theParticleIterator->value();
    G4ProcessManager* pManager = pParticle->GetProcessManager();
    G4String szParticleName = pParticle->GetParticleName();

    if (szParticleName == "gamma") {
    // gamma
      // Construct processes for gamma
      pManager->AddDiscreteProcess(new G4GammaConversion());
      pManager->AddDiscreteProcess(new G4ComptonScattering());      
      pManager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } else if (szParticleName == "e-") {
    //electron
      // Construct processes for electron
      pManager->AddProcess(new G4MultipleScattering(),-1,1,1);
      pManager->AddProcess(new G4eIonisation(),-1,2,2);
      pManager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);
  
    } else if (szParticleName == "e+") {
    //positron
      // Construct processes for positron
     pManager->AddProcess(new G4MultipleScattering(),-1,1,1);
     
     pManager->AddProcess(new G4eIonisation(),-1,2,2);
     pManager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);      
     pManager->AddProcess(new G4eplusAnnihilation(),0,-1,4);
  
    } else if( szParticleName == "mu+" || 
               szParticleName == "mu-"    ) {
    //muon  
     // Construct processes for muon+
     pManager->AddProcess(new G4MultipleScattering(),-1,1,1);
     pManager->AddProcess(new G4MuIonisation(),-1,2,2);
     pManager->AddProcess(new G4MuBremsstrahlung(),-1,-1,3);
     pManager->AddProcess(new G4MuPairProduction(),-1,-1,4);       
     
    }
    MyLDP* pNewLDP = new MyLDP;
    std::vector<G4VEvaporationChannel*>* MyChannels;
    MyChannels = new std::vector<G4VEvaporationChannel*>;
    MyChannels->reserve(8);
    ProtonEvaporationChannel* pProtonChannel = new ProtonEvaporationChannel();
    pProtonChannel->SetEmissionStrategy(new ProtonEvapProb());
    pProtonChannel->SetLevelDensityParameter(pNewLDP);
    MyChannels->push_back(pProtonChannel);
    NeutronEvaporationChannel* pNeutronChannel = new NeutronEvaporationChannel();
    pNeutronChannel->SetEmissionStrategy(new NeutronEvapProb);
    pNeutronChannel->SetLevelDensityParameter(pNewLDP);
    MyChannels->push_back(pNeutronChannel);
    DeuteronEvaporationChannel* pDeutChannel = new DeuteronEvaporationChannel();
    pDeutChannel->SetEmissionStrategy(new DeuteronEvapProb());
    pDeutChannel->SetLevelDensityParameter(pNewLDP);
    MyChannels->push_back(pDeutChannel); // Deuteron
    TritonEvaporationChannel* pTritonChannel = new TritonEvaporationChannel();
    pTritonChannel->SetEmissionStrategy(new TritonEvapProb());
    pTritonChannel->SetLevelDensityParameter(pNewLDP);
    MyChannels->push_back(pTritonChannel);   // Triton
    He3EvaporationChannel* pHe3Channel = new He3EvaporationChannel();
    pHe3Channel->SetEmissionStrategy(new He3EvapProb);
    pHe3Channel->SetLevelDensityParameter(pNewLDP);
    MyChannels->push_back(pHe3Channel);      // He3
    AlphaEvaporationChannel* pAlphaChannel = new AlphaEvaporationChannel();
    pAlphaChannel->SetEmissionStrategy(new AlphaEvapProb);
    pAlphaChannel->SetLevelDensityParameter(pNewLDP);
    MyChannels->push_back(pAlphaChannel);    // Alpha
    G4CompetitiveFission* pCompFission = new G4CompetitiveFission();
    pCompFission->SetLevelDensityParameter(new FissionLDP());
    MyChannels->push_back( pCompFission); // Fission Channel
    MyChannels->push_back( new G4PhotonEvaporation() );  // Photon Channel
    G4Evaporation* pMyEvaporation = new G4Evaporation(MyChannels);

    if(szParticleName=="proton"){
      G4hIonisation* pIonisation = new G4hIonisation();
      pManager->AddProcess(pIonisation);
      pManager->SetProcessOrdering(pIonisation,idxAlongStep);
      pManager->SetProcessOrdering(pIonisation,idxPostStep);
      pManager->AddProcess(new G4MultipleScattering(),-1,1,1);
      pManager->AddDiscreteProcess(pElasticProcess);
      pProtonLowEnergyModel = new G4ProtonInelasticProcess("ProtonLEinelastic");
      G4LEProtonInelastic* pProtonModel = new G4LEProtonInelastic;
      pProtonLowEnergyModel->RegisterMe(pProtonModel);


      pProtonPrecompoundModel = new G4ProtonInelasticProcess("ProtonPrecompound");
      G4ExcitationHandler* pProtonExcitation = new G4ExcitationHandler();
      G4PreCompoundModel* pProtonPrecompound = new G4PreCompoundModel(pProtonExcitation);
      pProtonPrecompoundModel->RegisterMe(pProtonPrecompound);

      //New stoi za Novata evaporacija. All - za syshtoto no za vsichki energii
      pProtonPrecompNew = new G4ProtonInelasticProcess("ProtonPrecompoundNew");
      pProtonPrecompAll = new G4ProtonInelasticProcess("ProtonPrecompoundAll");
      pProtonExcitation = new G4ExcitationHandler();
      pProtonExcitation->SetEvaporation(pMyEvaporation);
      pProtonExcitation->SetMinEForMultiFrag(1*GeV);
      pProtonPrecompNew->RegisterMe(new Model(pProtonExcitation));
      pProtonExcitation = new G4ExcitationHandler();
      pProtonExcitation->SetEvaporation(pMyEvaporation);
      pProtonExcitation->SetMinEForMultiFrag(1*GeV);
      pProtonPrecompAll->RegisterMe(new Model(pProtonExcitation));

      pProtonMars5GeVModel = new G4ProtonInelasticProcess("ProtonMars5GeV");
      G4Mars5GeV* pProtonMars = new G4Mars5GeV();
      pProtonMars->SetMinEnergy(0);
      pProtonMars5GeVModel->RegisterMe(pProtonMars);
      pProtonCutProcess = new Mars01EminCut();
      pProtonCutProcess->SetEminThreshold(1*MeV);

      pProtonChiralInvModel = new G4ProtonInelasticProcess("chiral invariant phase space for proton");
      G4PionMinusNuclearReaction* pProtChiral = new G4PionMinusNuclearReaction();
      pProtonChiralInvModel->RegisterMe(pProtChiral);

      G4CrossSectionDataStore* pStore = ((G4HadronInelasticProcess*)pProtonLowEnergyModel)->GetCrossSectionDataStore();
      G4CrossSectionDataStore* pStore_1 = ((G4HadronInelasticProcess*)pProtonPrecompoundModel)->GetCrossSectionDataStore();
      G4CrossSectionDataStore* pStore_2 = ((G4HadronInelasticProcess*)pProtonMars5GeVModel)->GetCrossSectionDataStore();
      G4ProtonInelasticCrossSection* protonData = new G4ProtonInelasticCrossSection;
      pStore->AddDataSet(protonData);
      pStore_1->AddDataSet(protonData);
      pStore_2->AddDataSet(protonData);
      ((G4HadronInelasticProcess*)pProtonPrecompNew)->GetCrossSectionDataStore()->AddDataSet(protonData);
      ((G4HadronInelasticProcess*)pProtonPrecompNew)->GetCrossSectionDataStore()->AddDataSet(protonData);

      pManager->AddDiscreteProcess(pProtonLowEnergyModel);
    }
    else if(szParticleName == "alpha"){
      pManager->AddProcess(new G4hIonisation(),-1,2,2);
      pManager->AddProcess(new G4MultipleScattering(),-1,1,1);
      
      pAlphaLowEnergy = new G4HadronInelasticProcess("alpha inelastic",pParticle);
      pAlphaLowEnergy->RegisterMe(new G4LEAlphaInelastic);

      pAlphaPrecompound = new G4HadronInelasticProcess("alpha precompound",pParticle);
      G4ExcitationHandler* pTypanar = new G4ExcitationHandler;
      pTypanar->SetMinEForMultiFrag(1*GeV);
      G4PreCompoundModel* pAlphaPrecomp = new G4PreCompoundModel(pTypanar);
      pAlphaPrecompound->RegisterMe(pAlphaPrecomp);

      pAlphaMars = new G4HadronInelasticProcess("alpha mars",pParticle);
      G4Mars5GeV* pAlphamars = new G4Mars5GeV;
      pAlphamars->SetMinEnergy(0);
      pAlphaMars->RegisterMe(pAlphamars);
      pAlphaCut = new Mars01EminCut();
      pAlphaCut->SetEminThreshold(1*MeV);

      pAlphaChiral = new G4HadronInelasticProcess("alpha chiral",pParticle);
      pAlphaChiral->RegisterMe(new G4PionMinusNuclearReaction);

      G4HadronElasticProcess* pAlphaElastic = new G4HadronElasticProcess;
      pAlphaElastic->RegisterMe(new G4LElastic);
      pManager->AddDiscreteProcess(pAlphaElastic);
      pManager->AddDiscreteProcess(pAlphaLowEnergy);
    }
    else if(szParticleName == "deuteron"){
      pManager->AddProcess(new G4MultipleScattering(),-1,1,1);
      pManager->AddProcess(new G4hIonisation(),-1,2,2);

      pDeuteronLowEnergy = new G4HadronInelasticProcess("deuteron inelastic",pParticle);
      pDeuteronLowEnergy->RegisterMe(new G4LEDeuteronInelastic);

      pDeuteronPrecompound = new G4HadronInelasticProcess("deuteron precompound",pParticle);
      G4ExcitationHandler* pTypanar2 = new G4ExcitationHandler;
      pTypanar2->SetMinEForMultiFrag(1*GeV);
      G4PreCompoundModel* pDeuteronPrecomp = new G4PreCompoundModel(pTypanar2);
      pDeuteronPrecompound->RegisterMe(pDeuteronPrecomp);

      pDeuteronMars = new G4HadronInelasticProcess("deuteron mars",pParticle);
      G4Mars5GeV* pDeutMars = new G4Mars5GeV;
      pDeutMars->SetMinEnergy(0);
      pDeuteronMars->RegisterMe(pDeutMars);
      pDeuteronCut = new Mars01EminCut;
      pDeuteronCut->SetEminThreshold(1*MeV);

      pDeuteronChiral = new G4HadronInelasticProcess("deuteron chiral",pParticle);
      pDeuteronChiral->RegisterMe(new G4PionMinusNuclearReaction);

      G4HadronElasticProcess* pElDeut = new G4HadronElasticProcess();
      pElDeut->RegisterMe(new G4LElastic);
      pManager->AddDiscreteProcess(pElDeut);
      pManager->AddDiscreteProcess(pDeuteronLowEnergy);
    }
    else if(szParticleName=="neutron"){
      G4LElastic* pElasticModel1 = new G4LElastic();
      G4NeutronHPElastic* pElasticNeutron = new G4NeutronHPElastic;
      pElasticProcess1->RegisterMe(pElasticModel1);
      pElasticModel1->SetMinEnergy(19.*MeV);
      pElasticProcess1->RegisterMe(pElasticNeutron);
      G4CrossSectionDataStore* pStore = ((G4HadronElasticProcess*)pElasticProcess1)->GetCrossSectionDataStore();
      G4NeutronHPElasticData* pNeutronData = new G4NeutronHPElasticData;
      pStore->AddDataSet(pNeutronData);
      pManager->AddDiscreteProcess(pElasticProcess1);

      pLowEnergyModel = new G4NeutronInelasticProcess("inelastic");
      G4LENeutronInelastic* pInelastic = new G4LENeutronInelastic;
      pInelastic->SetMinEnergy(19.*MeV);
      pLowEnergyModel->RegisterMe(pInelastic);
      m_cRegisteredModel = LEinelastic;

      pPrecompoundModel = new G4NeutronInelasticProcess("precompound");
      G4ExcitationHandler* pExcitation = new G4ExcitationHandler();
      G4PreCompoundModel* pPrecompound = new G4PreCompoundModel(pExcitation);
      pPrecompound->SetMinEnergy(19*MeV);
      pPrecompoundModel->RegisterMe(pPrecompound);

      //novija precompound
      pPrecompoundNew = new G4NeutronInelasticProcess("PrecompoundNew");
      pPrecompoundAll = new G4NeutronInelasticProcess("PrecompoundAll");
      pExcitation = new G4ExcitationHandler();
      pExcitation->SetEvaporation(pMyEvaporation);
      pExcitation->SetMinEForMultiFrag(1*GeV);
      Model* pPrecompoundM = new Model(pExcitation);
      pPrecompoundM->SetMinEnergy(19*MeV);
      pPrecompoundNew->RegisterMe(pPrecompoundM);
      pPrecompoundM = new Model(pExcitation);
      pPrecompoundM->SetMinEnergy(0*MeV);
      pPrecompoundAll->RegisterMe(pPrecompoundM);

      pMars5GeVModel = new G4NeutronInelasticProcess("Mars");
      G4Mars5GeV* pMarsProcess = new G4Mars5GeV();
      pMarsProcess->SetMinEnergy(19*MeV);
      pMars5GeVModel->RegisterMe(pMarsProcess);
      pCutProcess = new Mars01EminCut();
      pCutProcess->SetEminThreshold(1*MeV);

      pChiralInvModel = new G4NeutronInelasticProcess("Chiral invariant phase space for neutron");
      G4PionMinusNuclearReaction* pNeutChiral = new G4PionMinusNuclearReaction();
      pNeutChiral->SetMinEnergy(19*MeV);
      pChiralInvModel->RegisterMe(pNeutChiral);

      G4NeutronHPInelastic* pInelasticNeutronModel = new G4NeutronHPInelastic;
      pLowEnergyModel->RegisterMe(pInelasticNeutronModel);
      pPrecompoundModel->RegisterMe(pInelasticNeutronModel);
      pMars5GeVModel->RegisterMe(pInelasticNeutronModel);
      pChiralInvModel->RegisterMe(pInelasticNeutronModel);
      pPrecompoundNew->RegisterMe(pInelasticNeutronModel);
      G4CrossSectionDataStore* pStore1 = ((G4HadronInelasticProcess*)pLowEnergyModel)->GetCrossSectionDataStore();
      G4CrossSectionDataStore* pStore1_1 = ((G4HadronInelasticProcess*)pPrecompoundModel)->GetCrossSectionDataStore();
      G4CrossSectionDataStore* pStore1_2 = ((G4HadronInelasticProcess*)pMars5GeVModel)->GetCrossSectionDataStore();
      G4CrossSectionDataStore* pStore1_3 = ((G4HadronInelasticProcess*)pChiralInvModel)->GetCrossSectionDataStore();
      G4NeutronInelasticCrossSection * nDataPrecise = new G4NeutronInelasticCrossSection;
      pStore1->AddDataSet(nDataPrecise);
      pStore1_1->AddDataSet(nDataPrecise);
      pStore1_2->AddDataSet(nDataPrecise);
      pStore1_3->AddDataSet(nDataPrecise);
      G4NeutronHPInelasticData* pData1 = new G4NeutronHPInelasticData;
      pStore1->AddDataSet(pData1);
      pStore1_1->AddDataSet(pData1);
      pStore1_2->AddDataSet(pData1);
      pStore1_3->AddDataSet(pData1);
      ((G4HadronInelasticProcess*)pPrecompoundNew)->GetCrossSectionDataStore()->AddDataSet(nDataPrecise);
      ((G4HadronInelasticProcess*)pPrecompoundAll)->GetCrossSectionDataStore()->AddDataSet(nDataPrecise);
      pManager->AddDiscreteProcess(pLowEnergyModel);

      G4HadronFissionProcess* pFissionProcess = new G4HadronFissionProcess;
      G4LFission* pFissionModel = new G4LFission;
      pFissionModel->SetMinEnergy(20.*MeV);
      pFissionProcess->RegisterMe(pFissionModel);
      G4NeutronHPFission* pLEFissionModel = new G4NeutronHPFission;
      pFissionProcess->RegisterMe(pLEFissionModel);
      G4CrossSectionDataStore* pStore2 = ((G4HadronFissionProcess*)pFissionProcess)->GetCrossSectionDataStore();
      G4NeutronHPFissionData* pData2 = new G4NeutronHPFissionData;
      pStore2->AddDataSet(pData2);
      pManager->AddDiscreteProcess(pFissionProcess);

      G4HadronCaptureProcess* pCaptureProcess = new G4HadronCaptureProcess;
      G4LCapture* pCaptureModel = new G4LCapture;
      pCaptureModel->SetMinEnergy(20.*MeV);
      pCaptureProcess->RegisterMe(pCaptureModel);
      G4NeutronHPCapture* pLECapture = new G4NeutronHPCapture;
      pCaptureProcess->RegisterMe(pLECapture);
      G4CrossSectionDataStore* pStore3 = ((G4HadronCaptureProcess*)pCaptureProcess)->GetCrossSectionDataStore();
      G4NeutronHPCaptureData* pData3 = new G4NeutronHPCaptureData;
      pStore3->AddDataSet(pData3);
      pManager->AddDiscreteProcess(pCaptureProcess);
    }
  }
}

#include "G4Decay.hh"

void PhysicsList::ConstructGlobal()
{
  G4Decay* pDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while((*theParticleIterator)()){
    G4ParticleDefinition* pParticle = theParticleIterator->value();
    G4ProcessManager* pManager = pParticle->GetProcessManager();
    if(pDecayProcess->IsApplicable(*pParticle)){
      pManager->AddProcess(pDecayProcess);
      pManager->SetProcessOrdering(pDecayProcess,idxPostStep);
      pManager->SetProcessOrdering(pDecayProcess,idxAtRest);
    }
  }
}
#include "G4ParticleTable.hh"
void PhysicsList::SetCuts()
{
  SetCutsWithDefault();
}

void PhysicsList::ResetModel(G4ProcessManager* pProtManager,G4ProcessManager* pManager,
			     G4ProcessManager* pAlphaManager,G4ProcessManager* pDeuteronManager)
{
  G4ProcessVector* pProtVector = pProtManager->GetProcessList();
  G4ProcessVector* pNeutVector = pManager->GetProcessList();
  G4ProcessVector* pAlphaVector = pAlphaManager->GetProcessList();
  G4ProcessVector* pDeuteronVector = pDeuteronManager->GetProcessList();
  switch(m_cRegisteredModel){
  case LEinelastic:
    if(pNeutVector->contains(pLowEnergyModel))
      pManager->RemoveProcess(pLowEnergyModel);
    if(pProtVector->contains(pProtonLowEnergyModel))
      pProtManager->RemoveProcess(pProtonLowEnergyModel);
    if(pAlphaVector->contains(pAlphaLowEnergy))
      pAlphaManager->RemoveProcess(pAlphaLowEnergy);
    if(pDeuteronVector->contains(pDeuteronLowEnergy))
      pDeuteronManager->RemoveProcess(pDeuteronLowEnergy);
    break;
  case Precompound:
    if(pNeutVector->contains(pPrecompoundModel))
      pManager->RemoveProcess(pPrecompoundModel);
    if(pProtVector->contains(pProtonPrecompoundModel))
      pProtManager->RemoveProcess(pProtonPrecompoundModel);
    if(pAlphaVector->contains(pAlphaPrecompound))
      pAlphaManager->RemoveProcess(pAlphaPrecompound);
    if(pDeuteronVector->contains(pDeuteronPrecompound))
      pDeuteronManager->RemoveProcess(pDeuteronPrecompound);
    break;
  case PrecompoundNew:
    if(pNeutVector->contains(pPrecompoundNew))
      pManager->RemoveProcess(pPrecompoundNew);
    if(pProtVector->contains(pProtonPrecompNew))
      pProtManager->RemoveProcess(pProtonPrecompNew);
    if(pAlphaVector->contains(pAlphaPrecompound))
      pAlphaManager->RemoveProcess(pAlphaPrecompound);
    if(pDeuteronVector->contains(pDeuteronPrecompound))
      pDeuteronManager->RemoveProcess(pDeuteronPrecompound);
    break;
  case PrecompoundAll:
    if(pNeutVector->contains(pPrecompoundAll))
      pManager->RemoveProcess(pPrecompoundAll);
    if(pProtVector->contains(pProtonPrecompAll))
      pProtManager->RemoveProcess(pProtonPrecompAll);
    if(pAlphaVector->contains(pAlphaPrecompound))
      pAlphaManager->RemoveProcess(pAlphaPrecompound);
    if(pDeuteronVector->contains(pDeuteronPrecompound))
      pDeuteronManager->RemoveProcess(pDeuteronPrecompound);
    break;
  case Mars:
    if(pNeutVector->contains(pMars5GeVModel))
      pManager->RemoveProcess(pMars5GeVModel);
    if(pNeutVector->contains(pCutProcess))
      pManager->RemoveProcess(pCutProcess);
    if(pProtVector->contains(pProtonMars5GeVModel))
      pProtManager->RemoveProcess(pProtonMars5GeVModel);
    if(pProtVector->contains(pProtonCutProcess))
      pProtManager->RemoveProcess(pProtonCutProcess);
    if(pAlphaVector->contains(pAlphaMars))
      pAlphaManager->RemoveProcess(pAlphaMars);
    if(pAlphaVector->contains(pAlphaCut))
      pAlphaManager->RemoveProcess(pAlphaCut);
    if(pDeuteronVector->contains(pDeuteronMars))
      pDeuteronManager->RemoveProcess(pDeuteronMars);
    if(pDeuteronVector->contains(pDeuteronCut))
      pDeuteronManager->RemoveProcess(pDeuteronCut);
    break;
  case Chiral:
    if(pNeutVector->contains(pChiralInvModel))
      pManager->RemoveProcess(pChiralInvModel);
    if(pProtVector->contains(pProtonChiralInvModel))
      pProtManager->RemoveProcess(pProtonChiralInvModel);
    if(pAlphaVector->contains(pAlphaChiral))
      pAlphaManager->RemoveProcess(pAlphaChiral);
    if(pDeuteronVector->contains(pDeuteronChiral))
      pDeuteronManager->RemoveProcess(pDeuteronChiral);
  default: break;
  }
}

void PhysicsList::ChangeModel(char Which)
{
  if((pLowEnergyModel==NULL)||(pPrecompoundModel==NULL)){
    G4cout<<"Warning: ConstructProcess is not called"<<G4endl;
    return;
  }
  G4ParticleDefinition* pNeutr = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
  G4ParticleDefinition* pProt = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  G4ParticleDefinition* pAlpha = G4ParticleTable::GetParticleTable()->FindParticle("alpha");
  G4ParticleDefinition* pDeuteron = G4ParticleTable::GetParticleTable()->FindParticle("deuteron");
  if((pNeutr==NULL)||(pProt==NULL)){
    G4cout<<"Warning: It seems that neutron || proton is not yet declared"<<G4endl;
    return;
  }
  if((pAlpha==NULL)||(pDeuteron==NULL)){
    G4cout<<"Warning: It seems that alpha or deuteron is not yet declared"<<G4endl;
    return;
  }
  if(Which==m_cRegisteredModel) return;
  G4ProcessManager* pManager = pNeutr->GetProcessManager();
  G4ProcessManager* pProtManager = pProt->GetProcessManager();
  G4ProcessManager* pAlphaManager = pAlpha->GetProcessManager();
  G4ProcessManager* pDeuteronManager = pDeuteron->GetProcessManager();
  ResetModel(pProtManager,pManager,pAlphaManager,pDeuteronManager);
  switch(Which){
  case LEinelastic:
    pManager->AddDiscreteProcess(pLowEnergyModel);
    pProtManager->AddDiscreteProcess(pProtonLowEnergyModel);
    pAlphaManager->AddDiscreteProcess(pAlphaLowEnergy);
    pDeuteronManager->AddDiscreteProcess(pDeuteronLowEnergy);
    m_cRegisteredModel = LEinelastic;
    break;
  case Precompound:
    pManager->AddDiscreteProcess(pPrecompoundModel);
    pProtManager->AddDiscreteProcess(pProtonPrecompoundModel);
    pAlphaManager->AddDiscreteProcess(pAlphaPrecompound);
    pDeuteronManager->AddDiscreteProcess(pDeuteronPrecompound);
    m_cRegisteredModel = Precompound;
    break;
  case PrecompoundNew:
    pManager->AddDiscreteProcess(pPrecompoundNew);
    pProtManager->AddDiscreteProcess(pProtonPrecompNew);
    pAlphaManager->AddDiscreteProcess(pAlphaPrecompound);
    pDeuteronManager->AddDiscreteProcess(pDeuteronPrecompound);
    m_cRegisteredModel = PrecompoundNew;
    break;
  case PrecompoundAll:
    pManager->AddDiscreteProcess(pPrecompoundAll);
    pProtManager->AddDiscreteProcess(pProtonPrecompAll);
    pAlphaManager->AddDiscreteProcess(pAlphaPrecompound);
    pDeuteronManager->AddDiscreteProcess(pDeuteronPrecompound);
    m_cRegisteredModel = PrecompoundAll;
    break;
  case Mars:
    pManager->AddDiscreteProcess(pMars5GeVModel);
    pManager->AddDiscreteProcess(pCutProcess);
    pProtManager->AddDiscreteProcess(pProtonMars5GeVModel);
    pProtManager->AddDiscreteProcess(pProtonCutProcess);
    pAlphaManager->AddDiscreteProcess(pAlphaMars);
    pAlphaManager->AddDiscreteProcess(pAlphaCut);
    pDeuteronManager->AddDiscreteProcess(pDeuteronMars);
    pDeuteronManager->AddDiscreteProcess(pDeuteronCut);
    m_cRegisteredModel = Mars;
    break;
  case Chiral:
    pManager->AddDiscreteProcess(pChiralInvModel);
    pProtManager->AddDiscreteProcess(pProtonChiralInvModel);
    pAlphaManager->AddDiscreteProcess(pAlphaChiral);
    pDeuteronManager->AddDiscreteProcess(pDeuteronChiral);
    m_cRegisteredModel = Chiral;
    break;
  default:
    G4cout<<"Warning: Unknown process id: "<<(int)Which<<G4endl;
    return;
  }
  m_pGeom->ResetParallelGeometry();
}
