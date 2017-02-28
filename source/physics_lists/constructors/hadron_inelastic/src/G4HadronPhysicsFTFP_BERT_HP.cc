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
//------------------------------------------------------------------------
//
// Modified:
//
//------------------------------------------------------------------------
//
#include "G4HadronPhysicsFTFP_BERT_HP.hh"

#include <iomanip>   
#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4CrossSectionInelastic.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4LFission.hh"

#include "G4CrossSectionDataSetRegistry.hh"

#include "G4PhysListUtil.hh"
#include "G4Threading.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTFP_BERT_HP);

G4ThreadLocal G4HadronPhysicsFTFP_BERT_HP::ThreadPrivate* 
G4HadronPhysicsFTFP_BERT_HP::tpdata=0;

G4HadronPhysicsFTFP_BERT_HP::G4HadronPhysicsFTFP_BERT_HP(G4int)
    : G4VPhysicsConstructor("hInelastic FTFP_BERT_HP")
/*    , theNeutrons(0)
    , theBertiniNeutron(0)
    , theFTFPNeutron(0)
    , theHPNeutron(0)
    , thePion(0)
    , theBertiniPion(0)
    , theFTFPPion(0)
    , theKaon(0)
    , theBertiniKaon(0)
    , theFTFPKaon(0)
    , thePro(0)
    , theBertiniPro(0)
    , theFTFPPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)*/
    , QuasiElastic(false)
  /*, xsKaon(0)
    , xsNeutronCaptureXS(0)*/
{}

G4HadronPhysicsFTFP_BERT_HP::G4HadronPhysicsFTFP_BERT_HP(const G4String& name, G4bool quasiElastic)
    : G4VPhysicsConstructor(name) 
/*    , theNeutrons(0)
    , theBertiniNeutron(0)
    , theFTFPNeutron(0)
    , theHPNeutron(0)
    , thePion(0)
    , theBertiniPion(0)
    , theFTFPPion(0)
    , theKaon(0)
    , theBertiniKaon(0)
    , theFTFPKaon(0)
    , thePro(0)
    , theBertiniPro(0)
    , theFTFPPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)*/
    , QuasiElastic(quasiElastic)
  /*, xsKaon(0)
    , xsNeutronCaptureXS(0)*/
{}

void G4HadronPhysicsFTFP_BERT_HP::CreateModels()
{

  G4double minFTFP_pion = 3.0 * GeV;
  G4double maxBERT_pion = 12.0 * GeV;
  G4double minFTFP_kaon = 3.0 * GeV;
  G4double maxBERT_kaon = 12.0 * GeV;
  G4double minFTFP_proton = 3.0 * GeV;
  G4double maxBERT_proton = 12.0 * GeV;
  G4double minFTFP_neutron = 3.0 * GeV;
  G4double maxBERT_neutron = 12.0 * GeV;

  if(G4Threading::IsMasterThread()) {
    G4cout << G4endl
         << " FTFP_BERT_HP : new threshold between BERT and FTFP is over the interval " << G4endl
         << " for pions :   " << minFTFP_pion/GeV << " to " << maxBERT_pion/GeV  << " GeV" << G4endl
         << " for kaons :   " << minFTFP_kaon/GeV << " to " << maxBERT_kaon/GeV  << " GeV" << G4endl
         << " for proton :  " << minFTFP_proton/GeV << " to " << maxBERT_proton/GeV  << " GeV" << G4endl
         << " for neutron : " << minFTFP_neutron/GeV << " to " << maxBERT_neutron/GeV  << " GeV" << G4endl
         << G4endl;
  }

  tpdata->theNeutrons=new G4NeutronBuilder( true ); // Fission on
  tpdata->theFTFPNeutron=new G4FTFPNeutronBuilder(QuasiElastic);
  tpdata->theNeutrons->RegisterMe(tpdata->theFTFPNeutron);
  tpdata->theFTFPNeutron->SetMinEnergy(minFTFP_neutron);
  tpdata->theNeutrons->RegisterMe(tpdata->theBertiniNeutron=new G4BertiniNeutronBuilder);
  tpdata->theBertiniNeutron->SetMinEnergy(19.9*MeV);
  tpdata->theBertiniNeutron->SetMaxEnergy(maxBERT_neutron);
  tpdata->theNeutrons->RegisterMe(tpdata->theHPNeutron=new G4NeutronPHPBuilder);

  tpdata->thePro=new G4ProtonBuilder;
  tpdata->theFTFPPro=new G4FTFPProtonBuilder(QuasiElastic);
  tpdata->thePro->RegisterMe(tpdata->theFTFPPro);
  tpdata->theFTFPPro->SetMinEnergy(minFTFP_proton);
  tpdata->thePro->RegisterMe(tpdata->theBertiniPro=new G4BertiniProtonBuilder);
  tpdata->theBertiniPro->SetMaxEnergy(maxBERT_proton);

  tpdata->thePion=new G4PionBuilder;
  tpdata->theFTFPPion=new G4FTFPPionBuilder(QuasiElastic);
  tpdata->thePion->RegisterMe(tpdata->theFTFPPion);
  tpdata->theFTFPPion->SetMinEnergy(minFTFP_pion);
  tpdata->thePion->RegisterMe(tpdata->theBertiniPion=new G4BertiniPionBuilder);
  tpdata->theBertiniPion->SetMaxEnergy(maxBERT_pion);

  tpdata->theKaon=new G4KaonBuilder;
  tpdata->theFTFPKaon=new G4FTFPKaonBuilder(QuasiElastic);
  tpdata->theKaon->RegisterMe(tpdata->theFTFPKaon);
  tpdata->theFTFPKaon->SetMinEnergy(minFTFP_kaon);
  tpdata->theKaon->RegisterMe(tpdata->theBertiniKaon=new G4BertiniKaonBuilder);
  tpdata->theBertiniKaon->SetMaxEnergy(maxBERT_kaon);
 
  tpdata->theHyperon=new G4HyperonFTFPBuilder;
    
  tpdata->theAntiBaryon=new G4AntiBarionBuilder;
  tpdata->theAntiBaryon->RegisterMe(tpdata->theFTFPAntiBaryon=new G4FTFPAntiBarionBuilder(QuasiElastic));
}

G4HadronPhysicsFTFP_BERT_HP::~G4HadronPhysicsFTFP_BERT_HP()
{
  if (!tpdata) return;

  delete tpdata->theNeutrons;
  delete tpdata->theBertiniNeutron;
  delete tpdata->theFTFPNeutron;
  delete tpdata->theHPNeutron;

  delete tpdata->thePion;
  delete tpdata->theBertiniPion;
  delete tpdata->theFTFPPion;
    
  delete tpdata->theKaon;
  delete tpdata->theBertiniKaon;
  delete tpdata->theFTFPKaon;
    
  delete tpdata->thePro;
  delete tpdata->theBertiniPro;
  delete tpdata->theFTFPPro;    
   
  delete tpdata->theHyperon;
  delete tpdata->theAntiBaryon;
  delete tpdata->theFTFPAntiBaryon;

  delete tpdata; tpdata = 0;
} 

void G4HadronPhysicsFTFP_BERT_HP::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void G4HadronPhysicsFTFP_BERT_HP::ConstructProcess()
{
  if (tpdata == 0 ) tpdata = new ThreadPrivate;
  CreateModels();
  tpdata->theNeutrons->Build();
  tpdata->thePro->Build();
  tpdata->thePion->Build();
  tpdata->theKaon->Build();

  // --- Kaons ---
  tpdata->xsKaon = new G4ComponentGGHadronNucleusXsc();
  G4VCrossSectionDataSet * kaonxs = new G4CrossSectionInelastic(tpdata->xsKaon);
  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(kaonxs);
    
  tpdata->theHyperon->Build();
  tpdata->theAntiBaryon->Build();

  // --- Neutrons ---
  G4HadronicProcess* capture = 0;
  G4HadronicProcess* fission = 0;
  G4ProcessManager* pmanager = G4Neutron::Neutron()->GetProcessManager();
  G4ProcessVector*  pv = pmanager->GetProcessList();
  for ( size_t i=0; i < static_cast<size_t>(pv->size()); ++i ) {
    if ( fCapture == ((*pv)[i])->GetProcessSubType() ) {
      capture = static_cast<G4HadronicProcess*>((*pv)[i]);
    } else if ( fFission == ((*pv)[i])->GetProcessSubType() ) {
      fission = static_cast<G4HadronicProcess*>((*pv)[i]);
    }
  }
  if ( ! capture ) {
    capture = new G4HadronCaptureProcess("nCapture");
    pmanager->AddDiscreteProcess(capture);
  }
  tpdata->xsNeutronCaptureXS = (G4NeutronCaptureXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronCaptureXS::Default_Name());
  capture->AddDataSet(tpdata->xsNeutronCaptureXS);
  capture->AddDataSet( new G4ParticleHPCaptureData );
  G4NeutronRadCapture* theNeutronRadCapture = new G4NeutronRadCapture(); 
  theNeutronRadCapture->SetMinEnergy( 19.9*MeV ); 
  capture->RegisterMe( theNeutronRadCapture );
  if ( ! fission ) {
    fission = new G4HadronFissionProcess("nFission");
    pmanager->AddDiscreteProcess(fission);
  }
  G4LFission* theNeutronLEPFission = new G4LFission();
  theNeutronLEPFission->SetMinEnergy( 19.9*MeV );
  fission->RegisterMe( theNeutronLEPFission );
}
